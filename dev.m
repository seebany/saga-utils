%function [sigma_pwr, sigma_phi] =dev(rcvr_op,ph_truncated,pwr_truncated,sitenum_op,fluct)
% computes the standard deviation to be used for modeling Gaussian noise.
%
% Inputs: rcvr_op, an array of the receivers, just used for counting number
% of rxs.
% ph_truncated, cell array of the phase data for all receivers.
% pwr_truncated, cell array of the power data for all receivers. Output scale will be on amplitude, though, not power.
% sitenum_op
% fluct = 0 for phase or 1 for amplitude.
%
% Written by Aurora Lopez-Rubio 2019
% Commented by S. Datta-Barua 2 July 2020
function [sigma_pwr, sigma_phi] =dev(rcvr_op,ph_truncated,pwr_truncated,sitenum_op,fluct)
clear m
clear truncated_all
keyboard
%if fluct==0
    for rr = 1:size(rcvr_op, 1)
        if isempty(ph_truncated{rr})==0
            truncated_all(:,rr)=ph_truncated{rr};
        else
            scale=0.25;
            return
        end
    end

for i = 1:size(truncated_all,1)
    m(i)=mean(truncated_all(i,:));
end
A(rr+1,:)=m;
for rr = 1:size(rcvr_op,1)
    A(rr,:)=truncated_all(:,rr)';
    std_rr(rr,:) = std(A([rr,end],:));
end
for i = 1:size(std_rr,2)
    scale_rr(:,i) = mean(std_rr(:,i));
end
keyboard
sigma_phi = max(scale_rr) 
clear truncated_all

%elseif fluct==1      
    for rr = 1:size(rcvr_op, 1)
        if isempty(pwr_truncated{rr})==0
            truncated_all(:,rr)=pwr_truncated{rr};
        else
            scale=0.25;
            return
        end
    end
%end

for i = 1:size(truncated_all,1)
    m(i)=mean(truncated_all(i,:));
end
A(rr+1,:)=m;
for rr = 1:size(rcvr_op,1)
    A(rr,:)=truncated_all(:,rr)';
    std_rr(rr,:) = std(A([rr,end],:));
end
for i = 1:size(std_rr,2)
    scale_rr(:,i) = mean(std_rr(:,i));
end
sigma_pwr = max(scale_rr) 
% fig=figure
% rcvr_op2 = char(rcvr_op);
% rcvr_op2 = [rcvr_op2; 'gridmea'];
% for rr = 1:size(A,1)
%     plot(A(rr,:)', 'color', rx_color(rcvr_op2(rr, :)));
%     hold on;
% end
% legend([sitenum_op;'Mean signal'], 'orientation', 'horizontal');
% ylabel('Phase $\Phi_f$ [rad]');
% title('Shifted detrended Phase $\Phi_f$')
% tightfig;
% saveas(fig, '/data1/home/alopez35/mfigures/scale', 'png');
% close;

end
