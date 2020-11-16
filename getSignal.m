%function [signal, freq] = getSignal(signal_type)
% takes as input signal_type = 0, 2.
% and returns string signal, e.g., 'L1CA'
% and frequency freq in Hz.
% signal_types are numerical values as specified
% in CASES receiver documentation at
% http://cases.astraspace.net/documentation/scintdef.txt
% 
% Commented by S. Datta-Barua 1 July 2020

function [signal, freq] = getSignal(signal_type)
switch signal_type
    case 0
        signal = 'L1CA';
        freq = 1575.42 * 10^6;
        %         case 1
        %             signal = 'L2CM';
        %             freq = 1227.60 * 10 ^ 6;
    case 2
        signal = 'L2CL';
        freq = 1227.60 * 10^6;
        %         case 3
        %             signal = 'L2CLM';
        %             freq = 1227.60 * 10 ^ 6;
        %         case 4
        %             signal = 'L5I';
        %         case 5
        %             signal = 'L5Q';
        %         case 6
        %             signal = 'L5IQ';
        %         case 7
        %             signal = 'L1CA-ALT1';
        %             freq = 1575.42 * 10 ^ 6;
        %         case 8
        %             signal = 'CDMA-UHF-PILOT';
        %         case 9
        %             signal = 'CDMA-UHF-SYNC';
    otherwise
        signal = '';
        freq = NaN;
end
end
