dPhi plots distribution of deltaphi between leptons
ttbar plots ttbar system distributions


to create a new project:
./bin/ma5 --expert

to use the existing projects (dPhi and ttbar) run as follows:

cd dPhi/SampleAnalyzer
source setup.sh
make
make (sometimes have to make twice for some reason)

./SampleAnalyzer --analysis="dPhi" list.txt



These are the commands to get the powheg and MC@NLO lhe files (on lxplus):

cmsStage /store/lhe/2512/TTTo2l2Nu2B_7TeV_11000000events.lhe .

cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_00.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_01.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_02.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_03.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_04.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_05.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_06.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_07.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_08.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_09.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_10.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_11.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_12.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_13.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_14.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_15.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_16.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_17.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_18.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_19.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_20.lhe .
cmsStage /store/lhe/5217/ttbar_mcatnlo_7TeV_cteq6m_21.lhe .
