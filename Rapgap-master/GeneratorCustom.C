AliGenerator *
GeneratorCustom()
{


  AliGenExtFile *gener     = new AliGenExtFile(-1);
  AliGenReaderHepMC *reader = new AliGenReaderHepMC();

  // need to resolve files from grid, make sure we have a valid connection
    //if (!gGrid) TGrid::Connect("alien");
  
       // the first event number is calculated from the input here (startevent) and the number of events (CONFIG_NEVENTS)
        reader->SetFileName("test.hepmc"); // put the location of the input file
     
        gener->SetReader(reader); 
            
        return gener;
}
 
