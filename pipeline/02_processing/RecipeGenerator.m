%%
set(0,'defaultAxesFontSize',25)
clear; close all;
session_name="Postprocessed"

if ~isfolder(cd+"/"+session_name)
    mkdir(cd+"/"+session_name)
end

%%
%% Generate indexes

A={};B={};
A.s6=[1
2
3
4
5
6
7
8
9
10
11
12
13
14];
B.s6=[     1     2
     1     3
     1     7
     2     5
     2     8
     3     4
     3     8
     4     5
     4     9
     5     6
     5     9
     6     7
     6     9
     7     8];

A.s12=14+[1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
20
21
22
23
24
25
26
27
];
B.s12=9+[     1     3
     1     4
     1     5
     1     6
     1     8
     1     9
     2     3
     2     5
     2     6
     2     7
     2     8
     2     9
     3     5
     3     7
     3     8
     3     9
     4     5
     4     6
     4     7
     4     8
     4     9
     5     7
     5     9
     6     7
     6     8
     6     9
     7     8];
A.s24=41+[1
2
3
4
5
6
];
B.s24=18+[     1     2
     3     4
     5     6
     7     8
     9    10
    11    12];

C={};
C.CommunityIDX_Coal=[A.s6; A.s12; A.s24];
C.CommunityIDX_Sub1=[B.s6(:,1); B.s12(:,1); B.s24(:,1)];
C.CommunityIDX_Sub2=[B.s6(:,2); B.s12(:,2); B.s24(:,2)];



% TO do : similarity between replicates, TO do : relations between
% coalescence pairs, To do : relations between 

T= struct2table(C);
filename = 'Postprocessed/CoalescenceRecipe.xlsx';
writetable(T,filename,'Sheet',1);


A.N=[1
2
3
4
5
6
7
8
9
10
11
12
13
14
15];
B.N=[     1     2
     2     3
     3     4
     4     5
     5     6
     6     1
     1     4
     3     6
     1     3
     2     4
     3     5
     4     6
     5     1
     6     2
     2     5];

C={}
C.CommunityIDX_Coal_natural=[A.N];
C.CommunityIDX_Sub1_natural=[B.N(:,1)];
C.CommunityIDX_Sub2_natural=[B.N(:,2)];

T= struct2table(C);

writetable(T,filename,'Sheet',2);


%% Subcommunity




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
    Samplename="P"+string(fieldA')+"-"+arrayfun(@(x) sprintf('%02d',x), fieldB', 'UniformOutput', false);
    IDX=100*fieldA'+fieldB';
    Type={}
    Type(1:length(fieldA),1)={fni};
    CommunityIDX=(1:length(fieldA))';

    T_toadd=[array2table(Samplename), array2table(IDX), cell2table(Type), array2table(CommunityIDX)]
    T=[T; T_toadd]
end

%






