%%
set(0,'defaultAxesFontSize',25)
clear;
session_name="Postprocessed"

% addpath to any shared MATLAB helper functions, if needed.
if ~isfolder(cd+"/"+session_name)
    
    mkdir(cd+"/"+session_name)
end

%% Generate indexes

%type = "" % F(Final), T(Timeseries)
%communityType = % S(Synthetic), N(Natural)
%medium =  % L, M, H
%coalescenceType = % S(Sub) , C(Coalesced)
%communityIdx = % Any integer. For S-S : 1~9 is for 6 species, 10~18 is for 12 species, 19~30 is for 30 species
%replicate= %1,2

%type_medium_communityType_coalescenceIdx
%ex) F_S_L_S_1_1
A.F_S_L_S_1=[1*ones(1,9), 1*ones(1,9), 2*ones(1,9), 1*ones(1,3)];
B.F_S_L_S_1=[1:9, 24+(1:9), 48+(1:9), 48+(1:3)];
A.F_S_L_S_2=[1*ones(1,9), 1*ones(1,9), 2*ones(1,9), 1*ones(1,3)];
B.F_S_L_S_2=[12+(1:9), 36+(1:9), 60+(1:9), 51+(1:3)];

A.F_S_M_S_1=[8*ones(1,9), 8*ones(1,9), 2*ones(1,9), 1*ones(1,3)];
B.F_S_M_S_1=[48+(1:9), 72+(1:9), 72+(1:9), 54+(1:3)];
A.F_S_M_S_2=[8*ones(1,9), 8*ones(1,9), 2*ones(1,9), 1*ones(1,3)];
B.F_S_M_S_2=[60+(1:9), 84+(1:9), 84+(1:9), 57+(1:3)];

A.F_S_H_S_1=[2*ones(1,9), 2*ones(1,9), 3*ones(1,12)];
B.F_S_H_S_1=[(1:9), 24+(1:9), 72+(1:12)];
A.F_S_H_S_2=[2*ones(1,9), 2*ones(1,9), 3*ones(1,12)];
B.F_S_H_S_2=[12+(1:9), 36+(1:9), 84+(1:12)];

A.F_S_L_C_1=[4*ones(1,47)];
B.F_S_L_C_1=[(1:41) (43:48)];
A.F_S_L_C_2=[4*ones(1,47)];
B.F_S_L_C_2=97-[(1:41) (43:48)];

A.F_S_M_C_1=[5*ones(1,47)];
B.F_S_M_C_1=[(1:41) (43:48)];
A.F_S_M_C_2=[5*ones(1,47)];
B.F_S_M_C_2=97-[(1:41) (43:48)];


A.F_S_H_C_1=[6*ones(1,47)];
B.F_S_H_C_1=[(1:41) (43:48)];
A.F_S_H_C_2=[6*ones(1,47)];
B.F_S_H_C_2=97-[(1:41) (43:48)];

%Natural community
A.F_N_L_S_1=[3*ones(1,6)];
B.F_N_L_S_1=[1, 13, 25, 37, 49, 61];
A.F_N_L_S_2=[3*ones(1,6)];
B.F_N_L_S_2=1+[1, 13, 25, 37, 49, 61];

A.F_N_M_S_1=[3*ones(1,6)];
B.F_N_M_S_1=2+[1, 13, 25, 37, 49, 61];
A.F_N_M_S_2=[3*ones(1,6)];
B.F_N_M_S_2=3+[1, 13, 25, 37, 49, 61];

A.F_N_H_S_1=[3*ones(1,6)];
B.F_N_H_S_1=4+[1, 13, 25, 37, 49, 61];
A.F_N_H_S_2=[3*ones(1,6)];
B.F_N_H_S_2=5+[1, 13, 25, 37, 49, 61];

A.F_N_L_C_1=[7*ones(1,15)];
B.F_N_L_C_1=[1, 13, 25, 37, 49, 61, 73, 85 ,1+ [1, 13, 25, 37, 49, 61, 73]];
A.F_N_L_C_2=[7*ones(1,15)];
B.F_N_L_C_2=2+[1, 13, 25, 37, 49, 61, 73, 85 ,1+ [1, 13, 25, 37, 49, 61, 73]];

A.F_N_M_C_1=[7*ones(1,15)];
B.F_N_M_C_1=4+[1, 13, 25, 37, 49, 61, 73, 85 ,1+ [1, 13, 25, 37, 49, 61, 73]];
A.F_N_M_C_2=[7*ones(1,15)];
B.F_N_M_C_2=6+[1, 13, 25, 37, 49, 61, 73, 85 ,1+ [1, 13, 25, 37, 49, 61, 73]];

A.F_N_H_C_1=[7*ones(1,15)];
B.F_N_H_C_1=8+[1, 13, 25, 37, 49, 61, 73, 85 ,1+ [1, 13, 25, 37, 49, 61, 73]];
A.F_N_H_C_2=[7*ones(1,15)];
B.F_N_H_C_2=10+[1, 13, 25, 37, 49, 61, 73, 85 ,1+ [1, 13, 25, 37, 49, 61, 73]];


% Timeseries

% TO do : similarity between replicates, TO do : relations between
% coalescence pairs, To do : relations between 

T = table;
T= makeTable(T,A,B);
filename = 'Postprocessed/Metadata.xlsx';
save("Postprocessed/AB.mat", 'A', 'B')
writetable(T,filename,'Sheet',1);

%%




%%
function T=makeTable(T,A,B)
    fn = fieldnames(A);
    
    for i = 1:numel(fn)
        fni = string(fn(i));
        fieldA = A.(fni);
        fieldB = B.(fni);

        T=assignin(T, fni, fieldA,fieldB);
    end
end

function T=assignin(T, fni, fieldA, fieldB)
    SampleIDX="P"+string(fieldA')+"-"+arrayfun(@(x) sprintf('%02d',x), fieldB', 'UniformOutput', false);
    IDX=100*fieldA'+fieldB';
    Type={}
    Type(1:length(fieldA),1)={fni};
    CommunityIDX=string((1:length(fieldA))');


    Type_parsed=struct
    for i=1:length(fieldA)
        type=split(Type{i},'_');
        Type_parsed.Timepoint(i,1)=type(1);
        Type_parsed.CommunityOrigin(i,1)=type(2);
        Type_parsed.Medium(i,1)=type(3);
        Type_parsed.CoalescenceType(i,1)=type(4);
        Type_parsed.Replicate(i,1)=type(5);
    end

    T_toadd=[array2table(SampleIDX), array2table(IDX), struct2table(Type_parsed), array2table((CommunityIDX))]



    T=[T; T_toadd]
end

%






