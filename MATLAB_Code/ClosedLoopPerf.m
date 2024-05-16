function [stats, stimProps, params] = ClosedLoopPerf(fPath, varargin)
    % Analyzes the performance of the TMS-EEG closed-loop stimulation system
    % by examining the distribution of stimulation pulses across phase and
    % amplitude bins for each frequency band. The function also calculates
    % the Kullback-Leibler divergence between the phase/amplitude distributions
    % to identify the frequency band where stimulation pulses were delivered
    % with the greatest specificity for phase and amplitude.

    % Inputs:
    % fPath - path to the file containing the TMS-EEG data
    % Optional Inputs:
    % OFFSET - number of samples to offset the sampling of the EEG
    %          signal from the TMS pulse (default = 10) 
    % FREQRES - resolution of the frequency bands in octaves
    %           (default = 0.25)
    % PHASERES - resolution of the phase bins in degrees
    %            (default = 20)
    % BANDW - bandwidth of the filter frequency bands in octaves
    %         (default = 1.5)
    % AMPS - percentiles to use for amplitude binning
    %        (default = 0:20:100)
    % SAMPRATE - sampling rate of the EEG data in Hz
    %            (default = 1000)
    % PREWIN - time window before the TMS pulse to sample the 
    %          baseline EEG signal in seconds (default = 1)

    % Outputs:
    % stats - structure containing the following fields:
    %         KLVals - Kullback-Leibler divergence values for each 
    %                  frequency band
    %         BestFreq - frequency band with the greatest triggering 
    %                    specificity for phase/amplitude
    %         BestPhaseM - mean phase of the stimulation pulses in 
    %                      the best frequency band
    %         BestPhaseSD - standard deviation of the phase of the
    %                       stimulation pulses in the best frequency band
    %         BestPhaseR - circular correlation coefficient of the phase
    %                      of the stimulation pulses in the best frequency band
    %         BestAmp80 - percentage of stimulation pulses in the best
    %                     frequency band that were delivered at amplitudes
    %                     greater than the 80th percentile
    % stimProps - table containing the following columns:
    %             TimeIndex - index of the stimulation pulse in the EEG data
    %             Freq - frequency band of the stimulation pulse
    %             Epoch - epoch number of the stimulation pulse
    %             Amp - amplitude of the stimulation pulse
    %             Ph - phase of the stimulation pulse in degrees
    %             PhUw - unwrapped phase of the stimulation pulse relative
    %                    to the start of the stimulation epoch in degrees
    %             FreqInd - index of the frequency band
    %             AmpInd - index of the amplitude bin
    %             PhInd - index of the phase bin
    %             PhUwInd - index of the unwrapped phase bin
    % params - structure containing the following fields:
    %          FilePath - path to the file containing the TMS-EEG data
    %          Freqs - list of frequency bands
    %          Bandwidth - bandwidth of the frequency bands
    %          PhaseEdges - list of phase bin edges
    %          AmpEdges - list of amplitude bin edges
    %          FreqCents - list of central frequencies of the frequency bands

    % Example:
    % Process included example data file, sample 2 secods of EEG prior to
    % stimulation when estimating amplitude distribution.
    % [stats, stimProps, params] = ClosedLoopPerformance("FMT_TM_PAS1_001.lvm", ...
    %                               'PREWIN', 2)

    % process inputs
    if any(strcmp(varargin, 'OFFSET'))
        offset = varargin{find(strcmp(varargin,'OFFSET'))+1};
    else
        offset = 10;
    end

    if any(strcmp(varargin, 'FREQRES'))
        fRes = varargin{find(strcmp(varargin,'FREQRES'))+1};
    else
        fRes = 0.25;
    end

    if any(strcmp(varargin, 'PHASERES'))
        pRes = varargin{find(strcmp(varargin,'PHASERES'))+1};
    else
        pRes = 20;
    end

    if any(strcmp(varargin, 'BANDW'))
        bw = varargin{find(strcmp(varargin,'BANDW'))+1};
    else
        bw = 1.5;
    end

    if any(strcmp(varargin, 'AMPS'))
        ampBins = varargin{find(strcmp(varargin,'AMPS'))+1};
    else
        ampBins = 0:20:100;
    end

    if any(strcmp(varargin, 'SAMPRATE'))
        fs = varargin{find(strcmp(varargin,'SAMPRATE'))+1};
    else
        fs = 1000;
    end

    if any(strcmp(varargin,'PREWIN'))
        preWin = varargin{find(strcmp(varargin, 'PREWIN'))+1};
    else
        preWin = 1;
    end

    fList = 2.^(0:fRes:8);
    fWidth = bw/fRes;
    phEdges = -180:pRes:180;
    preSamps = fs*preWin;

    % load data
    data = readtable(fPath, 'FileType', 'text');
    stimEpochs = data.Var2;
    stimEpochs(stimEpochs~=1) = 0;
    stimTimes = data.Var3;
    stimTimes(stimTimes~=1) = 0;
    sig = data.Var4;
    
    % identify stimulation inds
    stimInds = find(diff(stimTimes)>0)-offset;
    
    % identify stimulation permissable epochs
    epochStarts = find(diff(stimEpochs)>0)-offset;
    epochStops = find(diff(stimEpochs)<0);
    
    if stimEpochs(1)
        epochStarts = [1; epochStarts];
    end
    if stimEpochs(end)
        epochStops = [epochStops; length(stimEpochs)];
    end
    
    if length(epochStarts) ~= length(epochStops)
        error('Stim epochs misaligned')
    else
        numEpochs = length(epochStarts);
    end
    
    stimProps = [];
    ampEdges = [];
    % process each frequency band
    for j = 1:(length(fList)-fWidth)
        centFreqs(j) = fList(j+round(fWidth/2));
        [b, a] = butter(2, fList([j j+fWidth])/(fs/2));
        sigFilt = filter(b,a,sig);
        sigHib = hilbert(sigFilt);
        sigAmp = abs(sigHib);
        sigPh = angle(sigHib);
    
        ampSamps = [];
        stimSamps = [];
        for k = 1:numEpochs
            indsEpoch = epochStarts(k):epochStops(k);
            sigAmpEpoch = sigAmp(indsEpoch);
            sigPhEpoch = sigPh(indsEpoch);
            sigPhEpochUw = unwrap(sigPhEpoch)/(2*pi);
            currStims = stimInds((stimInds>=epochStarts(k))&(stimInds<=epochStops(k)));
    
            if ~isempty(currStims)
                ampSamps = [ampSamps; sigAmp((currStims(1)-preSamps):currStims(1))];
                currStims = (currStims - epochStarts(k))+1;
                currStims = reshape(currStims,[numel(currStims),1]);
            else
                continue;
            end
            
            stimSamps = [stimSamps; ((currStims + epochStarts(k))), ...
                        centFreqs(j)*ones(size(currStims)), k*ones(size(currStims)), ...
                        sigAmpEpoch(currStims), 180*(sigPhEpoch(currStims)/pi), ...
                        sigPhEpochUw(currStims)];
        end

        ampEdges(:,j) = prctile(ampSamps,ampBins);
        stimSamps(:,end+1) = j;
        stimSamps(:,end+1) = discretize(stimSamps(:,4), ampEdges(:,j));
        stimSamps(:,end+1) = discretize(stimSamps(:, 5), phEdges);
        stimSamps(:,end+1) = round((stimSamps(:,6)+pi)/(2*pi));

        stimProps = [stimProps; stimSamps];
    end
    
    % calculate phase and amplitude histograms
    stimProps = num2cell(stimProps,1);
    stimProps = table(stimProps{:},'VariableNames', {'TimeIndex', 'Freq', ...
                      'Epoch', 'Amp', 'Ph', 'PhUw', ...
                      'FreqInd', 'AmpInd', 'PhInd', 'PhUwInd'});
    
    stimArray = table2array(stimProps(:,{'AmpInd', 'PhInd', 'FreqInd'}));
    stimHist = accumarray(stimArray, 1, [length(ampBins)-1, ...
                          length(phEdges)-1, length(centFreqs)]);
    
    phHist = squeeze(sum(stimHist, 1));
    ampHist = squeeze(sum(stimHist,2));


    % use KL Divergence to identify the frequency where stimulation pulses
    % were delivered with the greatest specificity for phase and amplitude
    for j = 1:length(centFreqs)
        klVals(j) = KLDiv(stimHist(:,:,j),[]);
    end
    [~, maxKLInd] = max(klVals);

    stimArray = table2array(stimProps(:,{'Amp', 'Ph', 'FreqInd'}));
    bestAmps = sort(stimArray(stimArray(:,3)==maxKLInd,1));
    bestPhases = pi*(stimArray(stimArray(:,3)==maxKLInd,2)/180);

    % calculate statistics
    stats.KLVals = klVals;
    stats.BestFreq = centFreqs(maxKLInd);
    stats.BestPhaseM = 180*circ_mean(bestPhases)/pi;
    stats.BestPhaseSD = 180*circ_std(bestPhases)/pi;
    stats.BestPhaseR = circ_r(bestPhases);
    stats.BestAmp80 = sum(bestAmps>ampEdges(4,maxKLInd))/length(bestAmps);

    % plot results
    fTicks = [1:10 20:10:100];
    phTicks = -180:45:180;
    figure('Name', 'Performance for ' + fPath, 'PaperOrientation','landscape', 'PaperPosition', [0.5 0.5 9 7]);

    % plot dependence of stimulation on phase/frequency
    subplot(2,4,1:3);
    [meshF, meshPh] = meshgrid(centFreqs, phEdges);
    phHist = [phHist; nan([1, size(phHist,2)])];
    pcH = pcolor(meshF, meshPh, phHist);
    set(pcH, 'LineStyle', 'none');
    set(gca, 'TickDir', 'out');
    set(gca, 'XScale', 'log');
    set(gca, 'XTick', fTicks);
    set(gca, 'Ytick', phTicks);
    line(centFreqs(maxKLInd)*[1 1], [-180 180], 'color', 'r')
    ylabel('Phase (degrees)');
    title('Phase distribution');
    grid on;
    set(gca,'Layer', 'top')

    subplot(2,4,4);
    phCentsRad = linspace(-180,180,size(phHist,1));
    barh(phCentsRad,phHist(:,maxKLInd),1,'r')
    hold on;
    plot(max(phHist(:,maxKLInd))*(cos(pi*phCentsRad/180)+1)/2,phCentsRad,'color','k', 'linestyle',':')
    set(gca, 'Ytick', phTicks);
    set(gca, 'TickDir', 'out');
    xlabel('Stim count');
    title('Phase tracking');
    
    subplot(2,4,5:7)
    [meshF, meshA] = meshgrid(centFreqs, ampBins);
    ampHist = [ampHist; nan([1,size(ampHist,2)])];
    pcH = pcolor(meshF, meshA, ampHist);
    set(pcH, 'LineStyle', 'none');
    set(gca, 'TickDir', 'out');
    set(gca, 'XScale', 'log');
    set(gca, 'XTick', fTicks);
    set(gca, 'Ytick', ampBins);
    line(centFreqs(maxKLInd)*[1 1], [-0 100], 'color', 'r')
    title('Amplitude distribution');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude (Percentile)');
    grid on;
    set(gca,'Layer', 'top')
    
    subplot(2,4,8);
    phCentsRad = linspace(-180,180,size(phHist,1));
    barh(ampBins(2:end),ampHist(1:(end-1),maxKLInd),1,'r')
    set(gca, 'Ytick', ampBins);
    set(gca, 'TickDir', 'out');
    xlabel('Stim count');
    title('Amplitude tracking');
    

    sampInds = (-500:500) + stimInds;
    meanWave = mean(sig(sampInds),1);
    figure('Name', 'Mean TMS waveform for ' + fPath); 
    plot((-500:500)/fs,meanWave)
    xlabel('Time relative to TMS (sec)')
    ylabel('Signal voltage (ADC units)')
    title('Mean TMS waveform')
    grid on;
    


    
    params = struct('FilePath', fPath, 'Freqs', fList, 'Bandwidth', fWidth, ...
                    'PhaseEdges', phEdges, 'AmpEdges', ampEdges, ...
                    'FreqCents', centFreqs);
end

function val = KLDiv(arr1,arr2)
% Calculates the Kullback-Leibler divergence between two arrays
% All values in arr1 and arr2 should be positive. Normalizes to
% probabilities (i.e. sum(arr1,'all')=1). If arr2 is empty then assigns it
% to uniform distribution version of arr1.

    if isempty(arr2)
        arr2 = mean(arr1,'all')*ones(size(arr1));
    end
   
    arr1 = arr1/sum(arr1,'all');
    arr2 = arr2/sum(arr2,'all');

    val = sum(arr1.*log(arr1./arr2), 'all', 'omitnan');
end