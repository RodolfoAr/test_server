
%addpath('\\som.surrey.ac.uk\Store\personal\pgr\RA00826\MATLAB\4.5.6-th04022019\dynare4.5.6t\matlab')

% RMSE


resume=1; % - 0 normal run - 1 load 
% Load pre-computed statistics  
% To be used in case the estimation procedure gets stuck.
% remember to adjust ZZ_start
ZZ_start=59; 
ZZ=100; % total number of samples
%ZZ_V=[ 5 7 10 11 13, 15:17 20 26 27, 31:34, 38:43, 45 47 55:56,];
FileSim=['RBCo2'];
FilterAlgo=['PF300'];
endo_var_names={'C' 'H' 'I' 'Z' 'mu' 'K' 'Y'};

modfile1       = 'RBCnoh_sim'; % DGP
modfile2       = 'RBCnoh_est'; % Estimation file
sim=0; % - 0 take pre-generated dataset - 1 generate dataset
est=1; % - 0 don't run estimation - 1 runs estimation

%-------------------------------------------------------------------------
ME=1; % ME=1 Appending measurement errors to the observables (needs est=1)
      % ME=2 Generate measurement errors vectors 
      %      (needs both sim=0 and est=0 not to edit previously simulated datasets
magnitude=0;

if magnitude==1 % small
    sigma_C0 = 0.004;
    sigma_H0 = 0.009;
    sigma_I0 = 0.0019;
else % large
    sigma_C0 = 0.004;
    sigma_H0 = 0.009;
    sigma_I0 = 0.0019;
end
%----------------------------------------------------
      
      
save_mode_check=1;  % number of mode_check windows


%--------------------------------------------------------
% Creating some useful string and objects
%--------------------------------------------------------
FileRunName=[FileSim,'_',FilterAlgo,'.mat'];

endo_est_end=size(endo_var_names,2); % real number of endogenous variables (used to handle size issues due to pruning)

T_ALL=NaN(ZZ,1);
RMSE_ALL=NaN(ZZ,8);
RMSE_obs=NaN(ZZ,1);
ARB_ALL=NaN(ZZ,8);
ARB_obs=NaN(ZZ,1);

ZZ_loop=ZZ+1;

%-------------------------------------------------------


    if sim==1
        for zz = ZZ_start:ZZ_loop

        %FileOriName=['temp_stoch_simul_ori_',num2str(zz)];    
        %FileObsName=['temp_stoch_simul_obs_',num2str(zz)];

        % Generate dataset

        %command_string = ['dynare ' modfile1 ' noclearall'];
        command_string = ['dynare ' modfile1];
        T = evalc(command_string);

        stoch_simul_ori_(ZZ_loop+1-zz).oo_.endo_simul=oo_.endo_simul;
        disp(zz)
        end
        
        save(FileSim,'stoch_simul_ori_');
    end 

    
% Generate ME errors
if ME==2
    load(FileSim);
for zz = ZZ_start:ZZ_loop
          stoch_simul_ori_(zz).oo_.C_ME=sigma_C0.*randn(size(stoch_simul_ori_(zz).oo_.endo_simul(1,:)',1),1);
          stoch_simul_ori_(zz).oo_.H_ME=sigma_H0.*randn(size(stoch_simul_ori_(zz).oo_.endo_simul(2,:)',1),1);
          stoch_simul_ori_(zz).oo_.I_ME=sigma_I0.*randn(size(stoch_simul_ori_(zz).oo_.endo_simul(3,:)',1),1);    
end 
        save(FileSim,'stoch_simul_ori_');
end
    
if est==1
      
      actual_=load(FileSim);
if resume==0
      est_results_(ZZ).oo_=[];
      est_results_(ZZ).M_=[];
      save('est_results_.mat', 'est_results_');
else
      load('est_results_.mat')
end

for zz = ZZ_start:ZZ_loop % ZZ_V %
    {zz datetime('now')}
      if ME==0
          C_obs=actual_.stoch_simul_ori_(zz).oo_.endo_simul(1,:)';
          H_obs=actual_.stoch_simul_ori_(zz).oo_.endo_simul(2,:)';
          I_obs=actual_.stoch_simul_ori_(zz).oo_.endo_simul(3,:)';
          save('temp_stoch_simul_obs','C_obs','H_obs','I_obs');
      elseif ME==1
          % additive measurement error
          C_obs=actual_.stoch_simul_ori_(zz).oo_.endo_simul(1,:)'+actual_.stoch_simul_ori_(zz).oo_.C_ME;
          H_obs=actual_.stoch_simul_ori_(zz).oo_.endo_simul(2,:)'+actual_.stoch_simul_ori_(zz).oo_.H_ME;
          I_obs=actual_.stoch_simul_ori_(zz).oo_.endo_simul(3,:)'+actual_.stoch_simul_ori_(zz).oo_.I_ME;
          save('temp_stoch_simul_obs','C_obs','H_obs','I_obs');
      end

            if zz>1
                tStart = tic;
            end

            % Resume output in case Dynare get stuck
            if resume==1
                load('resume.mat')
            end
            
            % Estimate
            command_string = ['dynare ' modfile2];
             T = evalc(command_string)
	         
             
            if zz>1 % eliminate the first simulation because the parallelization toolbox needs to be started

                % Efficiency
                tEnd = toc(tStart);
                T_ALL(zz-1)=tEnd;
                
                
                % ACCURACY - RMSE
                % all variables
                %actual_=load(FileOriName);
                endo_act=actual_.stoch_simul_ori_(zz).oo_.endo_simul(1:end,101:end); 
                endo_est=oo_.endo_simul(1:endo_est_end,101:end);

                RMSE_sample = (sqrt(mean((endo_act-endo_est).^2,2)));
                RMSE_sample_1 = sqrt(mean((endo_act(1,:)-endo_est(1,:)).^2));
                RMSE_sample_2 = sqrt(mean((endo_act(2,:)-endo_est(2,:)).^2));
                RMSE_sample_3 = sqrt(mean((endo_act(3,:)-endo_est(3,:)).^2));
                RMSE_sample_4 = sqrt(mean((endo_act(4,:)-endo_est(4,:)).^2));
                RMSE_sample_5 = sqrt(mean((endo_act(5,:)-endo_est(5,:)).^2));
                RMSE_sample_6 = sqrt(mean((endo_act(6,:)-endo_est(6,:)).^2));
                RMSE_sample_7 = sqrt(mean((endo_act(7,:)-endo_est(7,:)).^2));
                RMSE_sample_m = mean(sqrt(mean((endo_act-endo_est).^2,2)));
                
                %RMSE_mean=[RMSE_sample_1,RMSE_sample_2,RMSE_sample_3,RMSE_sample_4,RMSE_sample_5,RMSE_sample_6];
                %RMSE_sample_bis = mean(RMSE_mean);

                RMSE_ALL(zz-1,:)=[RMSE_sample_m RMSE_sample_1 RMSE_sample_2 RMSE_sample_3 RMSE_sample_4 RMSE_sample_5 RMSE_sample_6 RMSE_sample_7];


                % ACCURACY OF POINT ESTIMATE - ARB (Average Relative Bias)

                ARB_sample = sqrt(mean(((endo_act-endo_est)./endo_act).^2,2));
                ARB_sample_1 = sqrt(mean(((endo_act(1,:)-endo_est(1,:))./endo_act(1,:)).^2));
                ARB_sample_2 = sqrt(mean(((endo_act(2,:)-endo_est(2,:))./endo_act(2,:)).^2));
                ARB_sample_3 = sqrt(mean(((endo_act(3,:)-endo_est(3,:))./endo_act(3,:)).^2));
                ARB_sample_4 = sqrt(mean(((endo_act(4,:)-endo_est(4,:))./endo_act(4,:)).^2));
                ARB_sample_5 = sqrt(mean(((endo_act(5,:)-endo_est(5,:))./endo_act(5,:)).^2));
                ARB_sample_6 = sqrt(mean(((endo_act(6,:)-endo_est(6,:))./endo_act(6,:)).^2));
                ARB_sample_7 = sqrt(mean(((endo_act(7,:)-endo_est(7,:))./endo_act(7,:)).^2));
                ARB_sample_m = mean(sqrt(mean(((endo_act-endo_est)./endo_act).^2,2)));
                
                
                ARB_ALL(zz-1,:)=[ARB_sample_m ARB_sample_1 ARB_sample_2 ARB_sample_3 ARB_sample_4 ARB_sample_5 ARB_sample_6 ARB_sample_7];


                % only observables (NOT ADAPTED TO THIS CASE)
                endo_act1=actual_.stoch_simul_ori_(zz).oo_.endo_simul(4:end,101:end); 
                endo_est1=oo_.endo_simul(4:endo_est_end,101:end);

                RMSE_sample1 = sqrt(mean(sum((endo_act1-endo_est1).^2)));
                RMSE_obs(zz-1)=RMSE_sample1;

                ARB_sample1 = sqrt(mean(sum(((endo_act1-endo_est1)./endo_act1).^2)));
                ARB_obs(zz-1)=ARB_sample1;

                save('resume.mat', 'RMSE_ALL', 'ARB_ALL', 'RMSE_obs', 'ARB_obs', 'T_ALL');
                
                
                load('est_results_');
                est_results_(zz-1).oo_=oo_;
                est_results_(zz-1).M_=M_;
                est_results_(zz-1).mode=load('RBCnoh_est_mode.mat');
               
                save('est_results_.mat', 'est_results_');
                
                % save mode_check files
                zz_str=num2str(zz);
                mode_check=[FileSim '_' FilterAlgo '_' zz_str]
                
                if save_mode_check==2
                    %close all
                    mode_check=[FileSim '_' FilterAlgo '_' zz_str '.fig']
                    mode_check_b=[FileSim '_' FilterAlgo '_' zz_str '_b.fig']
                    copyfile('dsge_base2t_obs_CheckPlots1.fig', mode_check);
                    copyfile('dsge_base2t_obs_CheckPlots2.fig', mode_check_b);
                     
                  %saveas(gcf, mode_check, 'fig')

                else
                  saveas(gcf, mode_check, 'fig')
                end
                %close all
                
                
            end

 end
    
    {'T_ALL' 'RMSE_ALL' endo_var_names{:} 'RMSE_obs'}
    [T_ALL, RMSE_ALL, RMSE_obs]
    {'ARB_ALL' endo_var_names{:} 'ARB_obs'}
    [ARB_ALL  ARB_obs]

    T_FINAL=mean(T_ALL);

    fprintf('%d minutes and %f seconds\n', floor(T_FINAL/60), rem(T_FINAL,60));
    
    stats_sample_.T_ALL=T_ALL;
    stats_sample_.RMSE_ALL=RMSE_ALL;
    stats_sample_.ARB_ALL=ARB_ALL;
    
    stats_summary_.RMSE_FINAL=mean(RMSE_ALL)
    stats_summary_.ARB_FINAL=mean(ARB_ALL)
    stats_summary_.RMSE_FINAL_std=std(RMSE_ALL)
    stats_summary_.ARB_FINAL_max=max(ARB_ALL)

    stats_summary_.T_FINAL_min=floor(mean(T_ALL)/60);
    stats_summary_.T_FINAL_sec=rem(mean(T_ALL),60);
    stats_summary_.T_FINAL_tot=mean(T_ALL);
    
    save(FileRunName, 'est_results_','stats_summary_','stats_sample_')
   
end 

clear;


%simulated_data=[num2cell(oo_.endo_simul')];
%headers=[cellstr(M_.endo_names)', cellstr(M_.exo_names)'];
%xlswrite('trend_0_o1',[headers;simulated_data(1:2000,:)]);

%moments=[mean(oo_.endo_simul'); std(oo_.endo_simul'); skewness(oo_.endo_simul'); kurtosis(oo_.endo_simul')];


%load(temp_stoch_simul_ori)
%err = sqrt(mean(sum((X-EST).^2)));