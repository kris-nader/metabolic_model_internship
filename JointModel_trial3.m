ll = readCbModel(['iNF517.mat']);
sc = readCbModel(['iMM904.mat']);

%% restrict o2 for both
ll=changeRxnBounds(ll,'EX_o2_e',0,'l');
sc=changeRxnBounds(sc,'EX_o2_e',0,'l');
sc=setMedia(sc,'C:\Users\user\Desktop\saccharomyces_cerevisiae\CDM35','CDM35');
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
basalMedium = data(1:end,1);
compounds = data(1:end,2);
rxnsMedia = {};
for i = 1:length(basalMedium)
    rxnsMedia = union(rxnsMedia, strsplit(compounds{i},';'));
    rxnsPerCompound{i} = strcat('EX_',setdiff(strsplit(compounds{i},';'),''),'[u]');
end
rxnsMedia = strcat('EX_',setdiff(rxnsMedia,''),'[u]');

% 
% for all the elements in the JoingModel.rxns that are NOT part of the
% media and are for the joing model(so it starts with EX and not SCIEX or
% LLIEX
for i = 1:length(JointModel.rxns)
    if (ismember(JointModel.rxns(i),rxnsMedia)==0 & startsWith(JointModel.rxns(i),'EX_')==1)==1
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
    if ismember(exchange_aa(i),rxnsMedia)==0
        ind = getPosOfElementsInArray(exchange_aa(i),JointModel.rxns); 
        JointModel.lb(ind)=0;
    end
end


for i = 1:length(rxnsMedia)
    temp=getPosOfElementsInArray(rxnsMedia(i),JointModel.rxns);
    if isempty(temp)==0
        explicit_rxn=getRxn_cobraFormat(JointModel,temp);
        temp_str = split(explicit_rxn,"<=>");
        if strcmp(temp_str(2),' ')==1
            JointModel.lb(temp)=-10;
        else
            JointModel.ub(temp)=10;
        end  
    else 
        rxnsMedia(i);
        continue
    end  
end

%% make sure lactose is the carbon source

JointModel.lb(getPosOfElementsInArray({'EX_lcts[u]'},JointModel.rxns)) = -10;
JointModel.lb(getPosOfElementsInArray({'EX_glc__D[u]'},JointModel.rxns)) = 0;



% check oxygen supply for joint Model
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



% in replacement of YDCom.ub(YDCom.indCom.EXsp(:)) = 1000;
JointModel.ub(JointModel.indCom.EXsp(find(JointModel.indCom.EXsp(:,1)),1)) = 1000;%
JointModel.ub(JointModel.indCom.EXsp(find(JointModel.indCom.EXsp(:,2)),2)) = 1000;%


JointModel.lb(JointModel.indCom.EXcom)
JointModel.ub(JointModel.indCom.EXcom)

rxnBiomass = strcat(nameTagsModel, [ll.rxns(find(ll.c==1));sc.rxns(find(sc.c==1))]);  % biomass reaction names
rxnBiomassId = findRxnIDs(JointModel, rxnBiomass);  % ids
JointModel.infoCom.spBm = rxnBiomass;  % .spBm for organism biomass reactions
JointModel.indCom.spBm = rxnBiomassId;


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
printUptakeBoundCom(JointModel, 1);
