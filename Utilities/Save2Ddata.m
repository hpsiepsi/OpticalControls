function fName = Save2Ddata(settings, data, fspec)

   switch settings.expType
       case 105 % basic2DIR experiment
          s = int2str(floor(settings.cur_t2));
       otherwise 
          s = [int2str(floor(settings.cur_t2)) '_' int2str(floor(settings.expParam1))];
   end
   
   fName = [settings.outputbasepath '/' settings.molPrefix s settings.savePostfix];
   
   save(fName, 'data', 'fspec')

end