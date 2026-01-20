function [output] = asca_bkh(varargin)
%  [output] = asca_bkh(data, design, design_matrix, design_label, nperm, parent_dir, ...
% var_class, sam_class_leb,var_class_leb, design_for_plot);
% This function computes ANOVA, mafdr, and fdr and returns P values, FDR Corrected P, FDR, and Q values
%
% INPUTS:
% (1) data = data (double) e.g. X(Samples, Variables)
% (2) design = design for ANOVA e.g. [1 1 1 2 2 2 ... N N N]
% (3) design_matrix = e.g. [1 0; 0 1; 1 1] evaluate two factors and interaction between them
%
% OPTIONAL INPUTS:
% (4) design_label = e.g. {'treatment'} or   {'treatment' 'time'}
% (5) nperm  =  Number Permutations
% (6) parent_dir e.g. 'C\data'
% (7) var_class = ... ;
% (8) sam_class_leb = {{'fac1_lev1' 'fac1_lev2' 'fac1_lev3'}' {'fac2_lev1'
% 'fac2_lev2'}'};
% (9) var_class_leb = {{'fac1_lev1' 'fac1_lev2' 'fac1_lev3'}' {'fac2_lev1'
% 'fac2_lev2'}'};
% (10) design_for_plot = exactly like design, but if you also want 2-factor
% interaction term, provide this by including the 3rd class vector
% (e.g., temp-vs-time design)
%
% OUTPUTS:
% output - structure, which contains ASCA_results and generates plots
%
% 26.04.2021
% % %--------------------------
% % Bekzod Khakimov
% % Associate Professor
% % Chemometrics and Analytical Technology Research Group
% % Deprtment of Food Science, University of Copenhagen
% % Rolighedsvej 26, Frederiksberg, 1958, Denmark
% % Office: +45 3532-8184, Mobile: +45 2887-4454
% % Email: bzo@food.ku.dk

%% Get Inputs
data=varargin{1};
design=varargin{2};
design_matrix=varargin{3};
if nargin < 4 || isempty(varargin{4});  design_label='no_design';  else design_label=varargin{4}; end
if nargin < 5 || isempty(varargin{5});  nperm=2000;  else nperm=varargin{5}; end
if nargin < 6 || isempty(varargin{6});  parent_dir=cd;  else parent_dir=varargin{6}; end

if nargin<7 || isempty(varargin{7})
    var_class=ones(size(data,2),1);
else
    var_class=varargin{7};
end


if nargin<8 || isempty(varargin{8})
    for scl=1:size(design,2)
        u1=unique(design(:,scl));
        for i=1:length(u1)
            sam_class_leb2{i,1}=num2str(u1(i));
        end
        sam_class_leb{1,scl}=sam_class_leb2;
    end
else
    sam_class_leb=varargin{8};
end


if nargin<9 || isempty(varargin{9})
    u1=unique(var_class);
    for i=1:length(u1)
        var_class_leb{i,1}=num2str(u1(i));
    end
else
    var_class_leb=varargin{9};
end

if nargin<10 || isempty(varargin{10})
    design_for_plot=varargin{2};
else
    design_for_plot=varargin{10};
end



%% clean design_labels
% design_label(:,find(isletter(design_label)==0))=[];

%% Get Time
t=now;
d = datetime(t,'ConvertFrom','datenum')
time_now=strrep(char(d),':','-');

%% Create Directory
[s1 s2]=size(design_label);
if s2>1
    des2='';
    for i2=1:s2
        des2=[des2 '-' design_label{i2}];
        %des2=des2(:,1:end-1);
    end
    des2=des2(1,2:end);
else
    des2=design_label{1,1};
end
a3a=['\ASCA for ' des2 '-' time_now]; % change platform e.g., GCMS, AminAcid, Meta
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(parent_dir,a3a);
dirs_full=[parent_dir a3a];
cd(dirs_full);

%% Plot Design
figure,bar(design); axis tight;
xlabel('samples');
ylabel(des2);
lgd= legend(design_label,'Location','best','Box','off');
lgd.FontName = 'Times New Roman';
lgd.FontSize = 8;
title(['Study Design (data size is ' num2str(size(data,1)) ' x ' num2str(size(data,2)) ')']);
print([des2 'study-design'], '-dpng','-r0');
savefig([des2 'study-design']);
% close all

%% ASCA
data(isnan(data))=0;
[ s1] = scale_bkh(data,'auto');
s1(isnan(s1))=0;
 ASCA_results = asca_real_bkh(s1,design,design_label,nperm,design_matrix);
%  ASCA_results=ASCAcat(s1,design,design_label,nperm,design_matrix);
print([des2 '-ASCA-TABLE'], '-dpng','-r0');
savefig([des2 '-ASCA-TABLE']);
% close all

%
%
%
%% Plot ASCA-Scores
%
%
%

design_org=design_for_plot;

% cond2 = design_for_plot(:,2); %overlay1
% u2    = unique(cond2); %overlay1
% n2 = numel(u2); %overlay1
% col2  = lines(numel(u2)); %overlay1
% marker2 = {'o' 's' 'd' 'v' '^' '<' '>' 'p' 'h' 'o' 's'}; %overlay1
% marker2 = marker2(1:numel(u2)); %overlay1

for dd=1:size(design_for_plot,2)

    pi=dd;
    factor_x=dd;
    design=design_org(:,factor_x);

    u1=unique(design);
    marker1={'o' 's' 'd' 'v' '^' '<' '>' 'p' 'h' 'o' 's' 'd' 'v' '^' '<' '>' 'p' 'h' 'o' 's' 'd' 'v' '^' '<' '>' 'p' 'h'};
    marker1=repmat(marker1,1,10);
    col=color_bkh;
    pi=factor_x;
    figure
    asca_scores=ASCA_results.Effects{pi}.loads{1,1};
    asca_scores1{pi}=asca_scores;
    [c r]=size(asca_scores);
    if r~=1

        p2=design;
        p3=unique(p2);
        for pipi=1:length(p3)
            p4=find(p2==p3(pipi));
            scatter(asca_scores(p4,1), asca_scores(p4,2), 90, col(pipi,:), marker1{pipi}, 'filled');
            alpha_val = 0.7; % 0 = fully transparent, 1 = opaque
            set(gca, 'ColorOrder', col); % optional, if you use color cycling
            set(get(gca,'Children'), 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', alpha_val, 'MarkerEdgeAlpha', alpha_val);
            hold on;
            legend off
        end

        %overlay2

        %  h1 = gobjects(length(p3),1);
        % for pipi=1:length(p3)
        %     p4 = find(p2==p3(pipi));
        %     h1(pipi) = scatter( asca_scores(p4,1), asca_scores(p4,2), ...
        %                         90, col(pipi,:), marker1{pipi}, 'filled' );
        %     set(h1(pipi), 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor',[1 1 1]);
        %     hold on;
        % end

        %end of overlay2

        
        hold on;
        % Plot the Ellipse First
        for ki=1:length(u1)
            sc1=asca_scores(design==u1(ki),1);
            sc2=asca_scores(design==u1(ki),2);
            plot_ellips_bkh([sc1,sc2], 3, col(ki,:), 0.2);  hold on;
        end

          % overlay condition-2 as smaller, semi-transparent points
        % for k = 1:n2
        %     idx = (cond2 == u2(k));
        %     scatter( asca_scores(idx,1), ...
        %              asca_scores(idx,2), ...
        %              40, ...                 % smaller point‐size
        %              col2(k,:), ...
        %              marker2{k}, ...
        %              'filled', ...
        %              'MarkerFaceAlpha',0.5 );
        % end
        % 
        % h2 = gobjects(n2,1);
        % for k = 1:n2
        %     idx = (cond2 == u2(k));
        %     h2(k) = scatter( asca_scores(idx,1), asca_scores(idx,2), ...
        %                      40, col2(k,:), marker2{k}, 'filled', ...
        %                      'MarkerFaceAlpha', 0.5 );
        %     hold on;
        % end
        % 
        % % then add a separate legend for cond2
        % cond2_labels = cellfun(@(x)sprintf('C2=%s',x), sam_class_leb{1,2}, 'uni',0);
        % lgd2 = legend(cond2_labels,'Location','northeastoutside');
        % lgd2.FontSize = 7;
        %  % now make one legend combining cond-1 and cond-2 —
        % cond2_labels = cellfun(@(x) sprintf('C2=%s', x), sam_class_leb{1,2}, 'uni',0);
        % all_handles = [h1; h2];
        % all_labels  = [sam_class_leb{1,dd}; cond2_labels];
        % legend(all_handles, all_labels, ...
        %        'Location','northeastoutside', 'NumColumns', 2, 'FontSize', 7);

        %end of overlay3

         % --- 1) annotate each real point with its cond2 number ---
    % x = asca_scores(:,1);
    % y = asca_scores(:,2);
    % labels2 = arrayfun(@num2str, cond2, 'uni', false);
    % dx = (max(x)-min(x))*0.01;  dy = (max(y)-min(y))*0.01;
    % for i = 1:length(x)
    %     text(x(i)+dx, y(i)+dy, labels2{i}, ...
    %          'FontSize',8, 'HorizontalAlignment','left', ...
    %          'Color',[0 0 0]);
    % end
    % 
    % % --- 2) build dummy handles for cond2 legend entries ---
    % h2 = gobjects(n2,1);
    % for k = 1:n2
    %     h2(k) = scatter(NaN, NaN, ...
    %                     40, col2(k,:), marker2{k}, ...
    %                     'filled', 'MarkerFaceAlpha',0.5 );
    % end
    % % 
    % % --- 3) grab your cond1 labels & cond2 labels ---
    % c1_labels = sam_class_leb{1,dd};    % e.g. {'1';'2';'3'}
    % c2_labels = arrayfun(@(v) sprintf('C2=%s',num2str(v)), u2, 'uni', false);
    % 
    % % --- 4) one legend for both factors ---
    % allH    = [h1; h2];
    % allLabs = [c1_labels; c2_labels];
    % lgd = legend(allH, allLabs, ...
    %              'Location','northeastoutside', ...
    %              'NumColumns',2, ...
    %              'FontSize',7);

    %end of new overlay


        xlabel(['SC1 ' num2str(round(ASCA_results.Effects{pi}.varExp2(1)*100)) ' %']);
        ylabel(['SC2 ' num2str(round(ASCA_results.Effects{pi}.varExp2(2)*100)) ' %']);
        c1_labels = sam_class_leb{1,dd}; %overlay!!!!!!!
        lgd = legend(c1_labels, 'Location','best','Box','off'); %overlay!!!!!!!
        lgd.FontName = 'Times New Roman'; %overlay!!!!!!!
        lgd.FontSize = 8; %overlay!!!!!!!
        lgd=legend(sam_class_leb,'Location','best','BOX','OFF');
        lgd.FontName = 'Times New Roman';
        lgd.FontSize = 8;
        hold on;

    else

        p2=design;
        p3=unique(p2);
        for pipi=1:length(p3)
            p4=find(p2==p3(pipi));
            scatter(p4,asca_scores(p4,1), 90, col(pipi,:), marker1{pipi}, 'filled');
            alpha_val = 0.7; % 0 = fully transparent, 1 = opaque
            set(gca, 'ColorOrder', col); % optional, if you use color cycling
            set(get(gca,'Children'), 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', alpha_val, 'MarkerEdgeAlpha', alpha_val);
            hold on;
        end


        % Plot the Ellipse Second
        for pipi=1:length(p3)
            p4=find(p2==p3(pipi));
            plot_ellips_bkh([p4,asca_scores(p4,1)], 3, col(pipi,:), 0.2);
            hold on;
        end

        xlabel('samples');
        ylabel(['SC1 ' num2str(round(ASCA_results.Effects{pi}.varExp2(1)*100)) ' %']);
        legend(sam_class_leb{1,pi},'Location','best','BOX','OFF');
        hold on;
    end
    p_ef=[num2str(round((str2num(ASCA_results.ANOVAtab{factor_x+1,5})),6)) '/' num2str(round((str2num(ASCA_results.ANOVAtab{factor_x+1,3})*100),2)) '%)'];
   % title(['ASCA scores for ' design_label{1,factor_x} ' effect (p/ef = ' p_ef],'FontWeight','normal','FontSize',8);
    title([design_label{1,factor_x} ' effect (p/ef = ' p_ef],'FontWeight','normal','FontSize',8);
    fig_bkh
    print([design_label{1,factor_x} '-ASCA-Scores' ], '-dpng','-r0');
    savefig([design_label{1,factor_x} '-ASCA-Scores']);
    % close all


    %%% Plot ASCA-Loadings
    

    figure
    asca_loadings=ASCA_results.Effects{pi}.loads{1,2};
    asca_loadings1{pi}=asca_loadings;

    [c r]=size(asca_loadings);
    if r~=1
        var_class(isnan(var_class),:)=0;
        vu_l=length(unique(var_class));
        vu=unique(var_class);
        for vi=1:length(vu)
            v_id=find(var_class==vu(vi));
            if vu_l>1
                scatter(asca_loadings(v_id,1),asca_loadings(v_id,2), 50, col(vi+1,:), marker1{vi}, 'filled'); % CHANGE
            else
                scatter(asca_loadings(v_id,1),asca_loadings(v_id,2), 150, [0.35 0.35 0.35], 'o', 'filled');              
            end
            alpha_val = 1.0; % 0 = fully transparent, 1 = opaque
            set(gca, 'ColorOrder', col(vi,:)); % optional, if you use color cycling
            set(get(gca,'Children'), 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', alpha_val, 'MarkerEdgeAlpha', alpha_val);
            hold on;
        end
        xlabel(['SC1 ' num2str(round(ASCA_results.Effects{pi}.varExp2(1)*100)) ' %']);
        ylabel(['SC2 ' num2str(round(ASCA_results.Effects{pi}.varExp2(2)*100)) ' %']);
        lgd=legend(var_class_leb,'Location','best','BOX','OFF');
        lgd.FontName = 'Times New Roman';
        lgd.FontSize = 8;
        hold on;

        for ko=1:length(asca_loadings(:,1))
            text(asca_loadings(ko,1),asca_loadings(ko,2),num2str(ko),'FontSize',6); %CHANGE FONT size of text
        end

    else

        var_class(isnan(var_class),:)=0;
        vu_l=length(unique(var_class));
        vu=unique(var_class);
        for vi=1:length(vu)
            v_id=find(var_class==vu(vi));
            if vu_l>1
                scatter(v_id,asca_loadings(v_id,1), 150, col(vi,:), 'o', 'filled');
            else
                scatter(v_id,asca_loadings(v_id,1), 150, [0.35 0.35 0.35], 'o', 'filled');
            end
            alpha_val = 0.3; % 0 = fully transparent, 1 = opaque
            set(gca, 'ColorOrder', col); % optional, if you use color cycling
            set(get(gca,'Children'), 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', alpha_val, 'MarkerEdgeAlpha', alpha_val);
            hold on;
        end

        xlabel('variables');
        ylabel(['SC1 ' num2str(round(ASCA_results.Effects{pi}.varExp2(1)*100,2)) ' %']);
       % legend(var_class_leb,'Location','best','BOX','OFF');
        hold on;

        for ko=1:length(asca_loadings(:,1))
            text(ko,asca_loadings(ko,1),num2str(ko),'FontSize',6); hold on;
        end

    end
    p_ef=[num2str(round((str2num(ASCA_results.ANOVAtab{factor_x+1,5})),6)) '/' num2str(round((str2num(ASCA_results.ANOVAtab{factor_x+1,3})*100),2)) '%)'];
    % title(['ASCA loadings for ' design_label{1,factor_x} ' effect (p/ef = ' p_ef],'FontWeight','normal','FontSize',8);
    title([design_label{1,factor_x} ' effect (p/ef = ' p_ef],'FontWeight','normal','FontSize',8);
    fig_bkh
    print([design_label{1,factor_x} '-ASCA-Loadings'], '-dpng','-r0');
    savefig([design_label{1,factor_x} '-ASCA--Loadings']);
end

%% Save output
output.ASCA_results=ASCA_results;
output.data=data;
output.design=design_org;
output.design_label=design_label;
output.var_class=var_class;
save ASCA_results output

end

