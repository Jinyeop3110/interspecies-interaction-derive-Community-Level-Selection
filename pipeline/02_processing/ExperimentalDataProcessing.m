%% Load
set(0,'defaultAxesFontSize',22)

global MetadataGenerator

run('init_PH.m')
run('init_OD.m')
run('init_GC.m')
%run('init_CC_idx.m')

%% community time series


[OD.F_S_L_S_1, OD.F_S_L_S_2]=Merge_S_ODANDPH(LN_SC_OD);
[PH.F_S_L_S_1, PH.F_S_L_S_2]=Merge_S_ODANDPH(LN_SC_PH);
[OD.F_S_M_S_1, OD.F_S_M_S_2]=Merge_S_ODANDPH(MN_SC_OD);
[PH.F_S_M_S_1, PH.F_S_M_S_2]=Merge_S_ODANDPH(MN_SC_PH);
[OD.F_S_H_S_1, OD.F_S_H_S_2]=Merge_S_ODANDPH(HN_SC_OD);
[PH.F_S_H_S_1, PH.F_S_H_S_2]=Merge_S_ODANDPH(HN_SC_PH);

[OD.F_S_L_C_1, OD.F_S_L_C_2]=Merge_S_ODANDPH(LN_CC_OD);
[PH.F_S_L_C_1, PH.F_S_L_C_2]=Merge_S_ODANDPH(LN_CC_PH);
[OD.F_S_M_C_1, OD.F_S_M_C_2]=Merge_S_ODANDPH(MN_CC_OD);
[PH.F_S_M_C_1, PH.F_S_M_C_2]=Merge_S_ODANDPH(MN_CC_PH);
[OD.F_S_H_C_1, OD.F_S_H_C_2]=Merge_S_ODANDPH(HN_CC_OD);
[PH.F_S_H_C_1, PH.F_S_H_C_2]=Merge_S_ODANDPH(HN_CC_PH);



%%
[GC.F_S_L_S_1,GC.F_S_L_S_2]=Merge_S_GC(LN_SC_GC,Time);
[GC.F_S_M_S_1,GC.F_S_M_S_2]=Merge_S_GC(MN_SC_GC,Time);
[GC.F_S_H_S_1,GC.F_S_H_S_2]=Merge_S_GC(HN_SC_GC,Time);
[GC.F_S_L_C_1,GC.F_S_L_C_2]=Merge_S_GC(LN_CC_GC,Time);
[GC.F_S_M_C_1,GC.F_S_M_C_2]=Merge_S_GC(MN_CC_GC,Time);
[GC.F_S_H_C_1,GC.F_S_H_C_2]=Merge_S_GC(HN_CC_GC,Time);


%%

[OD.F_N_L_S_1, OD.F_N_L_S_2, OD.F_N_M_S_1, OD.F_N_M_S_2, OD.F_N_H_S_1, OD.F_N_H_S_2]=Merge_N_ODANDPH(NCS_SC_OD);
[PH.F_N_L_S_1, PH.F_N_L_S_2, PH.F_N_M_S_1, PH.F_N_M_S_2, PH.F_N_H_S_1, PH.F_N_H_S_2]=Merge_N_ODANDPH(NCS_SC_PH);
[GC.F_N_L_S_1, GC.F_N_L_S_2, GC.F_N_M_S_1, GC.F_N_M_S_2, GC.F_N_H_S_1, GC.F_N_H_S_2]=Merge_N_GC(NCS_SC_GC,Time);

[OD.F_N_L_C_1, OD.F_N_L_C_2,OD.F_N_M_C_1, OD.F_N_M_C_2,OD.F_N_H_C_1, OD.F_N_H_C_2]=Merge_N_ODANDPH(NCS_CC_OD);
[PH.F_N_L_C_1, PH.F_N_L_C_2, PH.F_N_M_C_1, PH.F_N_M_C_2, PH.F_N_H_C_1, PH.F_N_H_C_2]=Merge_N_ODANDPH(NCS_CC_PH);
[GC.F_N_L_C_1, GC.F_N_L_C_2, GC.F_N_M_C_1, GC.F_N_M_C_2, GC.F_N_H_C_1, GC.F_N_H_C_2]=Merge_N_GC(NCS_CC_GC,Time);


%%


load("Postprocessed/AB.mat")

T = table;
T= makeTable_withExperiment(T,A,B,OD,PH,GC);
filename = 'Postprocessed/Metadata.xlsx';
writetable(T,filename,'Sheet',1);

%%




%%
function T=makeTable_withExperiment(T,A,B,OD,PH,GC);

    fn = fieldnames(A);
    for i = 1:numel(fn)
        fni = string(fn(i));
        fieldA = A.(fni);
        fieldB = B.(fni);
        fieldOD = OD.(fni);
        fieldPH = PH.(fni);        
        fieldGC = GC.(fni);
        T=assignin(T, fni, fieldA,fieldB, fieldOD, fieldPH, fieldGC);
    end
end

function T=assignin(T, fni, fieldA, fieldB, fieldOD, fieldPH, fieldGC)
    

    SampleIDX="P"+string(fieldA')+"-"+arrayfun(@(x) sprintf('%02d',x), fieldB', 'UniformOutput', false);
    IDX=100*fieldA'+fieldB';
    Type={}
    Type(1:length(fieldA),1)={fni};
    CommunityIDX=string((1:length(fieldA))');


    Type_parsed=struct;
    for i=1:length(fieldA)
        type=split(Type{i},'_');
        Type_parsed.Timepoint(i,1)=type(1);
        Type_parsed.CommunityOrigin(i,1)=type(2);
        Type_parsed.Medium(i,1)=type(3);
        Type_parsed.CoalescenceType(i,1)=type(4);
        Type_parsed.Replicate(i,1)=type(5);
    end

    T_toadd=[array2table(SampleIDX), array2table(IDX), struct2table(Type_parsed), array2table((CommunityIDX)), array2table(fieldOD),array2table(fieldPH), array2table(fieldGC)]



    T=[T; T_toadd]
end


%%
function [O1,O2]=Merge_S_ODANDPH(S)
    A=S.s6;
    B=S.s12;
    C=S.s24;
    O1=[squeeze(A(:,1,:)); squeeze(B(:,1,:)) ; squeeze(C(:,1,:))];
    O2=[squeeze(A(:,2,:)); squeeze(B(:,2,:)) ; squeeze(C(:,2,:))];
end


function [O1,O2]=Merge_S_GC(S,Time)
    A=S.s6;
    B=S.s12;
    C=S.s24;
    O1_=[squeeze(A(:,1)); squeeze(B(:,1)) ; squeeze(C(:,1))];
    O2_=[squeeze(A(:,2)); squeeze(B(:,2)) ; squeeze(C(:,2))];
    
    for i=1:length(O1_)
        ODt=O1_{i};
        [O1(i,1) O1(i,2) O1(i,3)]=GRcalculate(ODt,Time);
    end
    for i=1:length(O2_)
        ODt=O2_{i};
        [O2(i,1) O2(i,2) O2(i,3)]=GRcalculate(ODt,Time);
    end
end


function [LN1, LN2, MN1, MN2, HN1, HN2]=Merge_N_ODANDPH(S)
    A=S.LN;
    B=S.MN;
    C=S.HN;
    LN1=squeeze(A(:,1,:));
    LN2=squeeze(A(:,2,:));
    MN1=squeeze(B(:,1,:)); 
    MN2=squeeze(B(:,2,:));
    HN1=squeeze(C(:,1,:)); 
    HN2=squeeze(C(:,2,:)); 
end


function [LN1, LN2, MN1, MN2, HN1, HN2]=Merge_N_GC(S,Time)
    A=S.LN;
    B=S.MN;
    C=S.HN;
    O1_=[squeeze(A(:,1)); squeeze(B(:,1)) ; squeeze(C(:,1))];
    O2_=[squeeze(A(:,2)); squeeze(B(:,2)) ; squeeze(C(:,2))];
    
    ODt_list=squeeze(A(:,1));
    for i=1:length(ODt_list)
        ODt=ODt_list{i};
        [LN1(i,1) LN1(i,2) LN1(i,3)]=GRcalculate(ODt,Time);
    end

    ODt_list=squeeze(A(:,2));
    for i=1:length(ODt_list)
        ODt=ODt_list{i};
        [LN2(i,1) LN2(i,2) LN2(i,3)]=GRcalculate(ODt,Time);
    end

    ODt_list=squeeze(B(:,1));
    for i=1:length(ODt_list)
        ODt=ODt_list{i};
        [MN1(i,1) MN1(i,2) MN1(i,3)]=GRcalculate(ODt,Time);
    end

    ODt_list=squeeze(B(:,2));
    for i=1:length(ODt_list)
        ODt=ODt_list{i};
        [MN2(i,1) MN2(i,2) MN2(i,3)]=GRcalculate(ODt,Time);
    end
        ODt_list=squeeze(C(:,1));
    for i=1:length(ODt_list)
        ODt=ODt_list{i};
        [HN1(i,1) HN1(i,2) HN1(i,3)]=GRcalculate(ODt,Time);
    end

    ODt_list=squeeze(C(:,2));
    for i=1:length(ODt_list)
        ODt=ODt_list{i};
        [HN2(i,1) HN2(i,2) HN2(i,3)]=GRcalculate(ODt,Time);
    end

end


%%
function output=FinalDay(input)
hold on
output=[]
temp=input.s6(:,:,end);
output=[output; temp(:)]
temp=input.s12(:,:,end);
output=[output; temp(:)]
temp=input.s24(:,:,end);
output=[output; temp(:)]
end


function Plot_pHvOD(data_pH,data_OD,color,dispname)
hold on

scatter(data_pH,data_OD,'filled','MarkerFaceColor',color,'DisplayName',dispname)
hold off
end

%%
function Plot_7daysTimeseries(data,color,dispname)
hold on
data=reshape(data,[],7);
for i=1:size(data,1)
    plot(1:7,data(i,:),'Color',color,'DisplayName',dispname);
end
hold off
end

function Plot_7daysTimeseries_FluctuationInference(data,f_list, color_list,dispname_list)
hold on
data=reshape(data,[],7);
for i=1:size(data,1)
    plot(1:7,data(i,:),'Color',color_list(f_list(i)+1,:),'DisplayName',dispname_list(f_list(i)+1));
end
hold off
end


function [p_OD_list, p_PH_list, f_list]=FluctuationInference(ODdata, PHdata, ODthreshold, PHthreshold)
ODdata=reshape(ODdata,[],7);

p_OD_list=arrayfun(@(x) isODFluctuation(ODdata(x,:)), 1:size(ODdata,1));
p_PH_list=arrayfun(@(x) isPHFluctuation(PHdata(x,:)), 1:size(ODdata,1));

f_list=((p_OD_list>ODthreshold)+(p_PH_list>PHthreshold))>0;
end

function p=isODFluctuation(dailyOD)
p=std(dailyOD(5:7))/mean(dailyOD(5:7));
end

function p=isPHFluctuation(dailyPH)
p=max(dailyPH(5:7))-min(dailyPH(5:7));
end


function output=Merge(input)
    output=[];
    temp=reshape(input.s6,[],7);
    output=[output; temp]
    temp=reshape(input.s12,[],7);
    output=[output; temp]
    temp=reshape(input.s24,[],7);
    output=[output; temp]
end


