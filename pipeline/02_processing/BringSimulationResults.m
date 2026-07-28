% Set to the directory holding the simulation .mat outputs before running.
data_path=getenv("COALESCENCE_SIMULATION");
%session_name='result_rand';
session_name='result_gaussian';
save_path='SimulationResult'+"/"+session_name;

if ~isfolder(save_path)
    mkdir(save_path)
end


Max_species_num=24*2;

SampleNum=100;
Snum_list=[6,12,24];
Interaction_list=[0.33,0.66,0.99];

T=table()

N_count=0;
for i=1:3
    for j=1:3
        N_count=N_count+1;
        T_toadd=struct

        fpath="Intstr_"+int2str(Snum_list(j))+"_hS_"+num2str(Interaction_list(i))+".mat"
        load(data_path+"/"+session_name+"/"+fpath)
        Snum=Snum_list(j)
        Interaction=Interaction_list(i)

        T_toadd.SampleIDX=("SM"+num2str(N_count)+"_")+num2str((1:100)','%.3d');
        T_toadd.CommunityIDX=(1:SampleNum)';

        T_toadd.Max_S= repmat(Max_species_num, [SampleNum,1]);
        T_toadd.S= repmat(Snum, [SampleNum,1]);
        T_toadd.I= repmat(Interaction, [SampleNum,1]);
        x1_list=[x1_list zeros([SampleNum,(Max_species_num-Snum*2)])];
        x1_list=x1_list./(sum(x1_list,2)+1e-10);
        x2_list=[x2_list zeros([SampleNum,(Max_species_num-Snum*2)])];
        x2_list=x2_list./(sum(x2_list,2)+1e-10);
        x3_list=[x3_list zeros([SampleNum,(Max_species_num-Snum*2)])];
        x3_list=x3_list./(sum(x3_list,2)+1e-10);

        T_toadd=[struct2table(T_toadd) array2table(x1_list) array2table(x2_list) array2table(x3_list)];
        T=[T; T_toadd];
    end
end


filename = save_path+'/processed_simulations.xlsx'
writetable(T,filename,'Sheet',1);

