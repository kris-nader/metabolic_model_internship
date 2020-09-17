ll = readCbModel(['iNF517.mat']);
sc = readCbModel(['iMM904.mat']);
sc=setMedia(sc,'C:\Users\user\Desktop\saccharomyces_cerevisiae\CDM35','CDM35');

changeCobraSolver('ibm_cplex', 'LP');
%% restrict o2 for both
ll=changeRxnBounds(ll,'EX_o2_e',0,'l');
sc=changeRxnBounds(sc,'EX_o2_e',0,'l');

%% find missing nutrients
FBAsolution_ll_0 = optimizeCbModel(ll,'max'); %0.0426
temp=findExcRxns(ll);
excRxns=ll.rxns(temp);
negative_flux=ll.rxns(find(FBAsolution_ll_0.x<0));
nutrients=intersect(excRxns,negative_flux);
positions=getPosOfElementsInArray(nutrients,ll.rxns);
lb_nutrients=ll.lb(positions);

ll=setMedia(ll,'C:\Users\user\Desktop\saccharomyces_cerevisiae\CDM35','CDM35');
FBAsolution_ll_1 = optimizeCbModel(ll,'max'); %0
allowedbackwards=ll.rxns(find(ll.lb<0));
temp=findExcRxns(ll);
excRxns=ll.rxns(temp);
allowedNutrients=intersect(allowedbackwards,excRxns);
missing=setdiff(nutrients,allowedNutrients);

%% restrict glucose
ll=changeRxnBounds(ll,'EX_glc__D_e',0,'l');

%% allow supplement of the media to be adenine,inosine,uracil
for i=[1,9,14]
   ll=changeRxnBounds(ll,missing{i},-10,'l');
end

%% apply regex to change format
sc.mets = regexprep(sc.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    sc.(fieldToBeCellStr{j})(cellfun(@isempty, sc.(fieldToBeCellStr{j}))) = {''};
end

ll.mets = regexprep(ll.mets, '_([^_]+)$', '\[$1\]');
% make all empty cells in cell arrays to be empty string
fieldToBeCellStr = {'metFormulas'; 'genes'; 'grRules'; 'metNames'; 'rxnNames'; 'subSystems'};
for j = 1:numel(fieldToBeCellStr)
    ll.(fieldToBeCellStr{j})(cellfun(@isempty, ll.(fieldToBeCellStr{j}))) = {''};
end

%% create the joint model
nameTagsModel = {'LL'; 'SC'};
JointModel = createMultipleSpeciesModel({ll; sc}, nameTagsModel);
JointModel.csense = char('E' * ones(1,numel(JointModel.mets)))';


%% make sure that the media is applied for the joint Model
[~, data] = xlsread('C:\Users\user\Desktop\saccharomyces_cerevisiae\CDM35','CDM35');
names = data(1:end,1);
compounds_fromMedia = data(1:end,2);
rxns_from_Media = {};
for i = 1:length(names)
    rxns_from_Media = union(rxns_from_Media, strsplit(compounds_fromMedia{i},';'));
    rxnsPerCompound{i} = strcat('EX_',setdiff(strsplit(compounds_fromMedia{i},';'),''),'[u]');
end
rxns_from_media1=rxns_from_Media
rxns_from_Media1 = strcat('EX_',setdiff(rxns_from_media1,''),'_e');
rxns_from_Media = strcat('EX_',setdiff(rxns_from_Media,''),'[u]');
% 
% for all the elements in the JointModel.rxns that are NOT part of the
% media and are for the joint model(so it starts with EX and not SCIEX or
% LLIEX
for i = 1:length(JointModel.rxns)
    if (ismember(JointModel.rxns(i),rxns_from_Media)==0 & startsWith(JointModel.rxns(i),'EX_')==1)==1
        ind = getPosOfElementsInArray(JointModel.rxns(i),JointModel.rxns); 
        JointModel.lb(ind) = -1;
    end
end

aa={'ala';'arg';'asn';'asp';'cys';'glu';'gln';'gly';'his';'ile';'leu';'met';'phe';'pro';
    'ser';'thr';'trp';'tyr';'val'};
exhange_aa={};
for i = 1:length(aa)
    exchange_aa(i) = strcat('EX_',setdiff(strsplit(aa{i},';'),''),'__L[u]');
end
exchange_aa(8)={'EX_gly[u]'};

for i =1:length(exchange_aa)
    if ismember(exchange_aa(i),rxns_from_Media)==0
        ind = getPosOfElementsInArray(exchange_aa(i),JointModel.rxns); 
        JointModel.lb(ind)=0;
    end
end

 % if its on the left we restrict to -10
for i = 1:length(rxns_from_Media)
    temp=getPosOfElementsInArray(rxns_from_Media(i),JointModel.rxns);
    if isempty(temp)==0
        explicit_rxn=getRxn_cobraFormat(JointModel,temp);
        temp_str = split(explicit_rxn,"<=>");
        if strcmp(temp_str(2),' ')==1
            JointModel.lb(temp)=-10;
        else
            JointModel.ub(temp)=10;
        end  
    else 
        rxns_from_Media(i);
        continue
    end  
end


%% make sure lactose is the carbon source
JointModel.lb(getPosOfElementsInArray({'EX_lcts[u]'},JointModel.rxns)) = -10;
JointModel.lb(getPosOfElementsInArray({'EX_glc__D[u]'},JointModel.rxns))=0;

%check oxygen supply for joint Model
temp=getPosOfElementsInArray({'EX_o2[u]'},JointModel.rxns);
explicit_rxn=getRxn_cobraFormat(JointModel,temp);
JointModel.lb(temp);
JointModel.lb(temp)=0;

temp=getPosOfElementsInArray({'SCIEX_o2[u]tr'},JointModel.rxns);
JointModel.lb(temp)=0;

temp=getPosOfElementsInArray({'LLIEX_o2[u]tr'},JointModel.rxns);
JointModel.lb(temp)=0;


% set nucleotides for joint model
temp=getPosOfElementsInArray({'EX_ade[u]'},JointModel.rxns);
explicit_rxn=getRxn_cobraFormat(JointModel,temp);
JointModel.lb(temp)=0;

temp=getPosOfElementsInArray({'EX_ins[u]'},JointModel.rxns);
explicit_rxn=getRxn_cobraFormat(JointModel,temp);
JointModel.lb(temp)=0;

temp=getPosOfElementsInArray({'EX_ura[u]'},JointModel.rxns);
explicit_rxn=getRxn_cobraFormat(JointModel,temp);
JointModel.lb(temp)=0;


[JointModel.infoCom, JointModel.indCom] = getMultiSpeciesModelId(JointModel, nameTagsModel);
disp(JointModel.infoCom);




rxns_from_Media1_com = regexprep(rxns_from_Media1,'_e$','[u]')
ind = getPosOfElementsInArray(rxns_from_Media1_com,JointModel.rxns);
JointModel.lb(JointModel.indCom.EXcom) = 0;
JointModel.lb(ind) = -10;





% in replacement of YDCom.ub(YDCom.indCom.EXsp(:)) = 1000;
JointModel.ub(JointModel.indCom.EXsp(find(JointModel.indCom.EXsp(:,1)),1)) = 1000;%
JointModel.ub(JointModel.indCom.EXsp(find(JointModel.indCom.EXsp(:,2)),2)) = 1000;%


JointModel.lb(JointModel.indCom.EXcom)
JointModel.ub(JointModel.indCom.EXcom)

rxnBiomass = strcat(nameTagsModel, [ll.rxns(find(ll.c==1));sc.rxns(find(sc.c==1))]);  % biomass reaction names
rxnBiomassId = findRxnIDs(JointModel, rxnBiomass);  % ids
JointModel.infoCom.spBm = rxnBiomass;  % .spBm for organism biomass reactions
JointModel.indCom.spBm = rxnBiomassId;

printUptakeBoundCom(JointModel, 1);


JointModel.c(rxnBiomassId) =1; 
fbaCom = optimizeCbModel(JointModel);
fbaCom.x(rxnBiomassId) 


options = struct();
options.GRguess = 13;  % initial guess for max. growth rate
options.GRtol = 1e-6;  % tolerance for final growth rate
options.algorithm = 1;  % use the default algorithm (simple guessing for bounds, followed by matlab fzero)

%to make it work for st, I changed line 56 in getCobraSolverParams: 
%Before, it was valDef.feasTol = 1e-9;. Now it is valDef.feasTol = 1e-6;
%basically, this parameter set an upper bound for the largest bound
%violation in the LP problem. 1e-6 is still good. 
% http://www-eio.upc.es/lceio/manuals/cplex-11/html/usrcplex/solveLP20.html

[sol, result] = SteadyCom(JointModel, options);

if ~strcmp(result.stat, 'infeasible')
    for jSp = 1:2
        fprintf('X_%s:  %.6f\n', JointModel.infoCom.spAbbr{jSp}, result.BM(jSp));
    end
    disp(result);
end

printUptakeBoundCom(JointModel, 1);
checkObjective(JointModel)
%% Graph

% percentage of maximum total biomass of the community required. 100 for sum(biomass) = 1 (1 is the default total biomass)
options.optBMpercent = 100;  
n = size(JointModel.S, 2);  % number of reactions in the model
% options.rxnNameList is the list of reactions subject to FVA. Can be reaction names or indices.
% Use n + j for the biomass variable of the j-th organism. Alternatively, use {'X_j'} 
% for biomass variable of the j-th organism or {'X_Ec1'} for Ec1 (the abbreviation in EcCom.infoCom.spAbbr)
options.rxnNameList = {'X_LL'; 'X_SC'};
options.optGRpercent = [89:0.2:99, 99.1:0.1:100];  % perform FVA at various percentages of the maximum growth rate, 89, 89.1, 89.2, ..., 100
[fvaComMin,fvaComMax] = SteadyComFVA(JointModel, options);

optGRpercentFBA = [89:2:99 99.1:0.1:100];  % less dense interval to save time because the results are always the same for < 99%
nGr = numel(optGRpercentFBA);
[fvaFBAMin, fvaFBAMax] = deal(zeros(numel(options.rxnNameList), nGr));
% change the objective function to the sum of all biomass reactions
JointModel.c(:) = 0;
JointModel.c(JointModel.indCom.spBm) = 1;
JointModel.csense = char('E' * ones(1, numel(JointModel.mets)));
s = optimizeCbModel(JointModel);  % run FBA
grFBA = s.f;
for jGr = 1:nGr
    fprintf('Growth rate %.4f :\n', grFBA * optGRpercentFBA(jGr)/100);
    [fvaFBAMin(:, jGr), fvaFBAMax(:, jGr)] = fluxVariability(JointModel, optGRpercentFBA(jGr), 'max', JointModel.infoCom.spBm, 2);
end

grComV = result.GRmax * options.optGRpercent / 100;  % vector of growth rates tested
lgLabel = {'{\itLL }';'{\itSC }'};
col = [235 135 255; 0 235 0 ]/255;  % color

% SteadyCom

hold on
x = [grComV(:); flipud(grComV(:))];
for j = 1:2
    y = [fvaComMin(j, :), fliplr(fvaComMax(j, :))];
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
    p(j, 1).Color = col(j, :);
    
end
tl(1) = title('\underline{SteadyCom}', 'Interpreter', 'latex');
% tl(1).Position = [0.7 1.01 0];
% ax(1) = gca;
% ax(1).XTick = 0.66:0.02:0.74;
% ax(1).YTick = 0:0.2:1;
% xlim([0.66 0.74])
% ylim([0 1])

lg = legend(lgLabel);
lg.Box = 'off';
yl(1) = ylabel('Relative abundance');
xl(1) = xlabel('Community growth rate (h^{-1})');




%% Analyze Pairwise Relationship Using SteadyComPOA
%options.savePOA = ['POA' filesep 'JointModel'];  % directory and fila name for saving POA results
options.optGRpercent = [99 90 70 50];  % analyze at these percentages of max. growth rate
% Nstep is the number of intermediate steps that the independent variable will take different values
% or directly the vector of values, e.g. Nsetp = [0, 0.5, 1] implies fixing the independent variable at the minimum,
% 50% from the min to the max and the maximum value respectively to find the attainable range of the dependent variable.
% Here use small step sizes when getting close to either ends of the flux range
a = 0.001*(1000.^((0:14)/14));
options.Nstep = sort([a (1-a)]);
[POAtable, fluxRange, Stat, GRvector] = SteadyComPOA(JointModel, options);


nSp = 2;
spLab = {'{\it LL }';'{\it SC }'};
mark = {'A'};
nPlot = 0;
for j = 1:nSp
    for k = 1:nSp
        if k > j
            nPlot = nPlot + 1;
            ax(j, k) = subplot(nSp-1, nSp-1, (k - 2) * (nSp - 1) + j);
            hold on
            for p = 1:size(POAtable{1, 1}, 3)
                x = [POAtable{j, j}(:, :, p);POAtable{j, j}(end:-1:1, :, p);...
                    POAtable{j, j}(1, 1, p)];
                y = [POAtable{j, k}(:, 1, p);POAtable{j, k}(end:-1:1, 2, p);...
                        POAtable{j, k}(1, 1, p)];
                plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2)
            end
            xlim([0.001 1])
            ylim([0.001 1])
            ax(j, k).XScale = 'log';
            ax(j, k).YScale = 'log';
            ax(j, k).XTick = [0.001 0.01 0.1 1];
            ax(j, k).YTick = [0.001 0.01 0.1 1];
            ax(j, k).YAxis.MinorTickValues=[];
            ax(j, k).XAxis.MinorTickValues=[];
            ax(j, k).TickLength = [0.03 0.01];
            xlabel(spLab{j});
            ylabel(spLab{k});
            tx(j, k) = text(10^(-5), 10^(0.1), mark{nPlot}, 'FontSize', 12, 'FontWeight', 'bold');
        end
    end
end
lg = legend(strcat(strtrim(cellstr(num2str(options.optGRpercent(:)))), '%'));
%lg.Position = [0.7246 0.6380 0.1700 0.2015];
lg.Box='off';
subplot(3, 3, 3, 'visible', 'off');
t = text(0.2, 0.8, {'% maximum';'growth rate'});
for j = 1:nSp
    for k = 1:nSp
        if k>j
            %ax(j, k).Position = [0.15 + (j - 1) * 0.3, 0.8 - (k - 2) * 0.3, 0.16, 0.17];
            %ax(j, k).Color = 'none';
        end
    end
end
