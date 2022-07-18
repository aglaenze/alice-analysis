AliAnalysisTaskNanoJPsi2016* AddNanoJPsi2016(Int_t period = 0, Bool_t MC = kFALSE)
{
  TString name = "NanoJPsi2016";
  // get the manager via the static access member. since it's static, you don't need
  // to create an instance of the class here to call the function
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method. 
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskNanoJPsi2016* task = new AliAnalysisTaskNanoJPsi2016(name.Data());   
    if(!task) return 0x0;
    // set some options
    task->SetPeriod(period);
    //mgr->GetMCtruthEventHandler()? task->SetMC(kTRUE) : task->SetMC(kFALSE);
	if (MC) task->SetMC(kTRUE);
	else task->SetMC(kFALSE);

    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
	mgr->ConnectOutput(task,1,mgr->CreateContainer("fOutputList", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("fOutputTree", TTree::Class(), AliAnalysisManager::kOutputContainer,fileName.Data()));
    //mgr->ConnectOutput(task,3,mgr->CreateContainer("fOutputSimTree", TTree::Class(), AliAnalysisManager::kOutputContainer,fileName.Data()));
	if (MC) mgr->ConnectOutput(task,3,mgr->CreateContainer("fOutputSimTree", TTree::Class(), AliAnalysisManager::kOutputContainer,fileName.Data()));
    
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
