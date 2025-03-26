function [ft_struct,ft_todss] = Robs_automaticICA(data,eye_channels_idx, NREMOVE)
%addpath([U '\NoiseTools\']);       % http://audition.ens.fr/adc/NoiseTools/

ft_struct = data ;

[B,A]=butter(2,1/(data.fsample/2), 'high');
x = cat(2,ft_struct.trial{:})';
%x = nt_trial2mat(eeg_asr.trial);
x(isnan(x)) = 0;
xoriginal = x;
tmp=nt_pca(filter(B,A,x(:,eye_channels_idx)));
mask=abs(tmp(:,1))>3*median(abs(tmp(:,1)));
C0=nt_cov(x);
C1=nt_cov(bsxfun(@times, x,mask));
[todss,pwr0,pwr1]=nt_dss0(C0,C1);

figure(98); clf; 
plot(pwr1./pwr0, '.-'); ylabel('score'); xlabel('component'); title ('eyeblink DSS');
eye_components=x*todss;
xfinal=nt_tsr(x,eye_components(:,1:NREMOVE));
figure(99); clf;
subplot 311; plot(eye_components(:,1:NREMOVE)); title('eye components');
subplot 312; plot(xfinal); title('eyeblinks removed');
subplot 313; plot(xoriginal); title('original data');

% Topographies
ft_dummy = struct;
ft_dummy.topo = todss;
for pc=1:size(todss,1)
    ft_dummy.label{pc,1} = ['PC' num2str(pc)];
end
ft_dummy.topolabel = ft_struct.label;
ft_todss = ft_dummy;

%%%% SEGMENT THE DATA BACK IN TRIALS %%%%
[nChan,nTimePoints] = size(ft_struct.trial{1,1});
ft_struct.trial = mat2cell(xfinal',nChan,cell2mat(cellfun(@(x) size(x,2),ft_struct.trial,'UniformOutput',0)));      

end
