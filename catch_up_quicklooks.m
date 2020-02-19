% Script to catch up on quicklook plotting.
% S. Datta-Barua
% 10 Apr 2019

cd ~/mfiles/saga/saga-utils/

for doy = 55:100
    
    cmdstr = ['main(2019,', num2str(doy) ', 1);'];
    eval(cmdstr)
end