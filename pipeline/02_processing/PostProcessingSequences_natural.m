%%
set(0,'defaultAxesFontSize',25)
clear; close all;
session_name="Postprocessed"

if ~isfolder(cd+"/"+session_name)
    mkdir(cd+"/"+session_name)
end

global Seqdir TAXAdir OTUdir


initialCutoff=200;
MAX_ASV_num=50;
A=readtable(Seqdir);

B=readtable(TAXAdir);
C=readtable(OTUdir);



A=A(:,1:initialCutoff+1);
[Samplename,Sequenced,OTUIDX]=Parse(A);

% For now, I manually remove OTU idx
OTUIDX_toremove=[];
[Sequenced_filtered,OTUIDX_filtered]=Filtering(Sequenced,OTUIDX,OTUIDX_toremove);
%

NormalizedAbundance=Sequenced_filtered./(sum(Sequenced_filtered,2)+1e-6);



E=[varfun(@(x) erase(x,"_F_filt.fastq.gz"),(Samplename)), array2table(NormalizedAbundance)];
E.Properties.VariableNames(1)="SampleIDX"
filename="Postprocessed/processed_Sequences_natural.xlsx";
writetable(E,filename,'Sheet',1)

B=B(OTUIDX_filtered,:);
C=C(OTUIDX_filtered,:);
D=[B,C(:,2)];
D.Properties.VariableNames(1)="ASV";
D.Properties.VariableNames(8)="UniqueSeuquences";
writetable(D,filename,'Sheet',2)


function T=assignin(T, fni, fieldA, fieldB)
    Samplename="P"+string(fieldA')+"-"+arrayfun(@(x) sprintf('%02d',x), fieldB', 'UniformOutput', false);
    IDX=100*fieldA'+fieldB';
    Type={}
    Type(1:length(fieldA),1)={fni};
    CommunityIDX=(1:length(fieldA))';

    T_toadd=[array2table(Samplename), array2table(IDX), cell2table(Type), array2table(CommunityIDX)]
    T=[T; T_toadd]
end



%% Generate indexes

function [Samplename,Sequenced, OTUIDX] = Parse(A)
    Samplename=A(:,1);
    Sequenced=A(:,2:end);
    Sequenced=table2array(Sequenced);
    OTUIDX=1:size(Sequenced,2);
end

function [Sequenced_filtered,OTUIDX_filtered] = Filtering(Sequenced,OTUIDX,OTUIDX_toremove)
    %IF rank of ASV > 80, or sum of detected ASV <500 then filterd out 
    %Filtering criteria : 0.5% for 

    Input_Sequences=sum(Sequenced(:));
    Sequenced(:,OTUIDX_toremove)=[];
    OTUIDX(OTUIDX_toremove)=[];


    OTU_cutoff_flag=(sum(Sequenced,1)<500) | (OTUIDX>200);
    Sequenced(:,OTU_cutoff_flag)=[];
    OTUIDX(OTU_cutoff_flag)=[];

    TT=Sequenced./sum(Sequenced,1);
    Shannon_metric=-sum(log(TT+1e-20).*TT);
    OTU_cutoff_flag=Shannon_metric<log(2);
    OTUIDX(OTU_cutoff_flag)=[];

    zeroflag=(Sequenced./(sum(Sequenced,2)+1e-6))<0.001;
    Sequenced(zeroflag)=0;
    Sequenced_filtered=Sequenced;
    OTUIDX_filtered=OTUIDX;
    Final_Sequences=sum(Sequenced_filtered(:));
    disp("we attained : "+sprintf("%f",Final_Sequences/Input_Sequences) + " of reads")
end

function [NormalizedAbundance, PresentFlag] = Processing(Sequenced_filtered,OTUIDX_filtered)
    NormalizedAbundance=Sequenced_filtered./(sum(Sequenced_filtered,2)+1e-6);
    PresentFlag=(NormalizedAbundance>0);
end
%
