#' DATA DICTIONARY
#' creates  print class and levels and whatnot
#' @param dataset A data frame name.
#' @param NAcode A quoted string ie "NA" or whatever your NA code is.
#' @param printalllevels TRUE OR FALSE print all levels of your variables
#' @return names and column number of variables or matches.
#' @examples
#' data(Bladder_Cancer)
#' dictionary(Bladder_Cancer, "NA", printalllevels=FALSE)
#'
#' @export

dictionary<-function(dataset,NAcode=NA,printalllevels=FALSE){
  dataset <- sjlabelled::remove_all_labels(dataset)
  classvec<-namevec<-numlevel<-numnull<-levelsvec<-rep(NA,dim(dataset)[2]);
  for(i in 1:dim(dataset)[2]) {
    
    classvec[i]<-class(dataset[,i])[1];
    namevec[i]<-names(dataset)[i];
    tempnumlevel<-"";
    if(class(dataset[,i])[1]=="factor"){
      numlevel[i]<-length(levels(dataset[,i]));
      for(j in 1:numlevel[i]){
        templevel<-as.character(levels(dataset[,i])[j]);
        if(printalllevels==F){tempnumlevel<-c(tempnumlevel,templevel);}
        if(printalllevels==T){tempnumlevel<-paste(tempnumlevel,templevel);}
        tempnumlevel<-tempnumlevel
      }
      
      if(printalllevels==F){ levelsvec[i]<-paste(tempnumlevel[2],", ",tempnumlevel[3],", ",tempnumlevel[4],", ",tempnumlevel[5],", ",tempnumlevel[6],sep="");}
      if(printalllevels==T){ levelsvec[i]<-tempnumlevel;}
      
      
      numnull[i]<-sum(I(dataset[,i]%in%c(NAcode,"")));
    }#end if FACTOR
    
    j<-1; tempnumlevel<-"";
    if(class(dataset[,i])[1]=="character"){  #print("we have a character!");
      numlevel[i]<-length(unique(dataset[,i]));
      for(j in 1:numlevel[i]){
        templevel<-as.character(unique(dataset[,i])[j]);
        tempnumlevel<-c(tempnumlevel,templevel);
      }
      #levelsvec[i]<-tempnumlevel;
      tempnumlevel<- unique(tempnumlevel);
      #print(tempnumlevel)
      if(printalllevels==F){ levelsvec[i]<-paste(tempnumlevel[2],", ",tempnumlevel[3],", ",tempnumlevel[4],", ",tempnumlevel[5],", ",tempnumlevel[6],sep="");}
      if(printalllevels==T){ levelsvec[i]<-tempnumlevel;}
      
      numnull[i]<-sum(I(dataset[,i]%in%c(NAcode,"",NA)));
    }
    
    if(class(dataset[,i])[1]=="logical"){  #print("we have a Logical!");
      numlevel[i]<-2;#length(unique(dataset[,i]));
      levelsvec[i]<-"True or False";
      numnull[i]<-sum(I(dataset[,i]%in%c(NAcode,"",NA)));
    }
    
    if(class(dataset[,i])[1]%in%c("numeric","integer")){  #print("we have a Logical!");
      xbar<-round(mean(dataset[,i],na.rm=T),3);
      temprange<-range(dataset[,i],na.rm=T);
      numlevel[i]<-"-";
      levelsvec[i]<-paste("mean = ",xbar,
                          " Range: ","(",temprange[1],", ",temprange[2],")",sep="");#length(unique(dataset[,i]));
      #numnull[i]<-sum(I(dataset[,i]%in%c(NAcode,"",NA)));
      numnull[i]<-sum(is.na(dataset[,i]));
    }
    
    if( (class(dataset[,i])[1]%in%c("Date","POSIXct","POSIXt"))) {
      #print("WHAT's THE PROBLEM")
      #print(names(dataset)[i])
      levelsvec[i]<-paste( " Range: ","(", min(dataset[,i], na.rm = TRUE)  ,", ", max(dataset[,i], na.rm = TRUE),")",sep="");#length(unique(dataset[,i]));
      #numnull[i]<-table(is.na(dataset[,i]))[2];
      numnull[i]<-sum(is.na(dataset[,i]));
      numlevel[i]<-"-";
      
    }
    if(!(class(dataset[,i])[1]%in%c("character","factor","logical","Date","POSIXct","POSIXt","numeric","integer"))) {
      print("WHAT's THE PROBLEM")
      print(names(dataset)[i])
      #print(summary(dataset[,i])[4])
      xbar<-summary(dataset[,i])[4]#round(mean(dataset[,i],na.rm=T),2);
      # print(xbar)#
      # temprange<-range(dataset[,i],na.rm=T);
      # print(temprange)
      numlevel[i]<-"-";
      #print(numlevel)
      levelsvec[i]<-paste("mean = ",xbar,
                          " Range: ","(",summary(dataset[,i])[1],", ",summary(dataset[,i])[6],")",sep="");#length(unique(dataset[,i]));
      numnull[i]<-table(is.na(dataset[,i]))[2];
    }
  }
  tempdf<-data.frame(namevec,classvec,numnull,numlevel,levelsvec)
  names(tempdf)<-c("Var","Class","NumMissing","NumLevels","Levels")
  return(tempdf);}
