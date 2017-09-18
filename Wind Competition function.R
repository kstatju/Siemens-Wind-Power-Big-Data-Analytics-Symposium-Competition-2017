

##################################
## Code for creating variable represents Error occurring "Before", "On", "After" Visit

Code_B_A <- function(x){
  a = as.POSIXct(x[1], format = '%Y-%m-%d %H:%M:%S')
  b = as.POSIXct(x[2], format = '%Y-%m-%d %H:%M:%S')

  if (a > b){
    return('1Before')
  }else if ((a <= b) & ((a + as.numeric(x[3]) * 60) >= b)){
    return('2In')
  }else if ((a + as.numeric(x[3]) * 60) < b){
    return('3After')
  }
}

#################################


############################################
## concate EventWarningStop, ManualStop and StopUrgency Variables

concatvar <- function(x, septerator = '.'){
  if(!require("stringr", quietly =T))  install.packages("stringr")
  require(stringr)
  return(paste0(str_trim(x[!is.na(x) & !(x %in% "")]), collapse = septerator))
  
}


#############################################


########################################################
## Code for Markov chain matrix

mcmat <- function(x, timevar, resvar, visit = 'VisitId',byvar = NULL){
  
  
  if (is.null(byvar)){
    
    dis_co = unlist(unique(x[resvar]))
    uni_visit = unlist(unique(x[visit]))
    len1 = length(uni_visit)
    len = length(dis_co)
    mcmatx = array(0, c(len, len), dimnames = list(dis_co, dis_co))
    for (j in 1:len1){
      x1 = x[x[visit] == uni_visit[j],]
      x1 = x1[order(x1[timevar], decreasing = FALSE),]
      n = dim(x1)
      if (n[1]>1){
        for (i in 2:n[1]){
          mcmatx[as.character(x1[resvar][i-1,]), as.character(x1[resvar][i,])] = mcmatx[as.character(x1[resvar][i-1,]), as.character(x1[resvar][i,])] + 1
        }
      }
    }
  }else{
    dis_co = unlist(unique(x[resvar]))
    len = length(dis_co)
    bydis_co = unlist(unique(x[byvar]))
    len1 = length(bydis_co)
    mcmatx = array(0, c(len, len, len1), dimnames = list(dis_co, dis_co, bydis_co))
    for (j in 1:len1){
      x1 = x[x[byvar]== as.character(bydis_co[j]),]
      uni_visit = dis_co = unlist(unique(x1[visit]))
      len2 = length(uni_visit)
      for (k in 1:len2){
        x2 = x1[x1[visit]== uni_visit[k],]
        x2 = x2[order(x2[timevar], decreasing = FALSE),]
        n = dim(x2)
        if (n[1]>1){
          for (i in 2:n[1]){
            mcmatx[as.character(x2[resvar][i-1,]), 
                   as.character(x2[resvar][i,]), 
                   as.character(x2[byvar][i,])] = mcmatx[as.character(x2[resvar][i-1,]), 
                                                         as.character(x2[resvar][i,]), 
                                                         as.character(x2[byvar][i,])] + 1
          }
        }
      }
    }
  }
  return(mcmatx)
  
}


##############################################################################
# find common patter based on Markoc Chain

patBasedOnMC <- function(x1, timevar, resvar, firstn = 5, inout = 2, visit = 'VisitId', byvar = 'Code_B_A'){
  if(!require("reshape2", quietly =T))  install.packages("reshape2")
  require(reshape2)
  require(dplyr)
    uni_visit = unlist(unique(x1[visit]))
    uni_BA = unlist(unique(x1[byvar]))
    
    new_frm = data.frame(matrix(vector(), 0, 17), stringsAsFactors = FALSE)
    names(new_frm) = c('Park_Name', 'FactorA', 'FactorB', 'FactorC', 'FactorD', 'StationID', 'VisitType', 
                       'VisitId', 'VisitDurMinutes' ,"bpattern", 'ipattern', 'apattern', 'pattern_all_n',
                       'common_pattern_all', 'common_pattern_n', 'common_BA_pattern_n', 'common_BA_pattern_all')
    l = 0
    for (i in 1:length(uni_visit)){
      l= l+1
      x2 = x1[x1[visit]== uni_visit[i],]
      mc = mcmat(x2, timevar = timevar, resvar = resvar, visit = visit)
      mc = mc/sum(mc)
      mc1 = setNames(melt(mc), c('Source', 'Target', 'Weight'))
      mc1 = mc1[mc1$Weight != 0,]
      mc1 = mc1[order(mc1$Weight, decreasing = T),]
      uni_cod = paste0(sort(unique(unlist(mc1[1:firstn, inout]))), collapse = ", " ,sep = '')
      
      comm_p = NULL
      comm_p_n = NULL
      biapatt = NULL
      biapattba_n = NULL
      biapattba = NULL
      for (j in 1:length(uni_BA)){
        x3 = x2[x2[byvar]== uni_BA[j],]
        if(dim(x3)[1]>0){
          mc = mcmat(x3, timevar = timevar, resvar = resvar, visit = visit)
          mc = mc/sum(mc)
          mc2 = setNames(melt(mc), c('Source', 'Target', 'Weight'))
          mc2 = mc2[mc2$Weight != 0,]
          mc2 = mc2[order(mc2$Weight, decreasing = T),]
          buni_cod = paste0(sort(unique(unlist(mc2[1:firstn, inout]))), collapse = ", " ,sep = '')
          biapatt = append(biapatt, buni_cod)
          if (j > 1){
            comm_p_n = intersect(comm_p_n, unique(unlist(mc2[1:firstn, inout])))
            comm_p = as.integer(intersect(comm_p, unique(unlist(mc2[,inout]))))
            if (uni_BA[j] == '1Before' || uni_BA[j] == '2In'){
              biapattba_n = intersect(biapattba_n, unique(unlist(mc2[1:firstn, 1])))
              biapattba = as.integer(intersect(biapattba, unique(unlist(mc2[,1]))))
            }
          }else {
            comm_p_n = unique(unlist(mc2[1:firstn, inout]))
            comm_p = unique(unlist(mc2[,inout]))
            biapattba_n = unique(unlist(mc2[1:firstn,1]))
            biapattba = unique(unlist(mc2[,1]))
          }
        }else{
          biapatt = append(biapatt, NA)
        }
      }
      comm_p_n = paste0(comm_p_n, collapse = ', ')
      comm_p = paste0(sort(comm_p), collapse = ', ')
      biapattba_n = paste0(sort(biapattba_n), collapse = ', ')
      biapattba = paste0(sort(biapattba), collapse = ', ')
      new_frm[l,] = c(as.character(x2$Park_Name[1]), as.character(x2$FactorA[1]), as.character(x2$FactorB[1]),
                      as.character(x2$FactorC[1]), as.character(x2$FactorD[1]), x2$StationID[1], 
                      as.character(x2$VisitType[1]), x2$VisitId[1], x2$VisitDurMinutes[1], biapatt, uni_cod,
                      comm_p, comm_p_n, biapattba_n, biapattba)
    }
    return(new_frm)
}


##############################################################
# Common pattern based on count


patBasedOnFreq <- function(x1, resvar, firstn = 5, visit = 'VisitId', byvar = 'Code_B_A'){
  if(!require("reshape2", quietly =T))  install.packages("reshape2")
  require(reshape2)
  require(dplyr)
  uni_visit = unlist(unique(x1[visit]))
  uni_BA = unlist(unique(x1[byvar]))
  
  new_frm = data.frame(matrix(vector(), 0, 17), stringsAsFactors = FALSE)
  names(new_frm) = c('Park_Name', 'FactorA', 'FactorB', 'FactorC', 'FactorD', 'StationID', 'VisitType', 
                     'VisitId', 'VisitDurMinutes' ,"bpattern", 'ipattern', 'apattern', 'pattern_all_n',
                     'common_pattern_all', 'common_pattern_n', 'common_BA_pattern_n', 'common_BA_pattern_all')
  l = 0
  for (i in 1:length(uni_visit)){
    l= l+1
    x2 = x1[x1[visit]== uni_visit[i],]
    mc = as.data.frame(table(x2[resvar]))
    mc$Freq = mc$Freq/sum(mc$Freq)
    mc1 = mc
    mc1 = mc1[mc1$Freq != 0,]
    mc1 = mc1[order(mc1$Freq, decreasing = T),]
    uni_cod = paste0(sort(unique(unlist(mc1[1:firstn, 1]))), collapse = ", " ,sep = '')
    
    comm_p = NULL
    comm_p_n = NULL
    biapatt = NULL
    biapattba_n = NULL
    biapattba = NULL
    for (j in 1:length(uni_BA)){
      x3 = x2[x2[byvar]== uni_BA[j],]
      if(dim(x3)[1]>0){
        mc = as.data.frame(table(x3[resvar]))
        mc$Freq = mc$Freq/sum(mc$Freq)
        mc2 = mc
        mc2 = mc2[mc2$Freq != 0,]
        mc2 = mc2[order(mc2$Freq, decreasing = T),]
        buni_cod = paste0(sort(unique(unlist(mc2[1:firstn, 1]))), collapse = ", " ,sep = '')
        biapatt = append(biapatt, buni_cod)
        if (j > 1){
          comm_p_n = intersect(comm_p_n, unique(unlist(mc2[1:firstn, 1])))
          comm_p = as.integer(intersect(comm_p, unique(unlist(mc2[,1]))))
          if (uni_BA[j] == '1Before' || uni_BA[j] == '2In'){
            biapattba_n = intersect(biapattba_n, unique(unlist(mc2[1:firstn, 1])))
            biapattba = as.integer(intersect(biapattba, unique(unlist(mc2[,1]))))
          }
        }else {
          comm_p_n = unique(unlist(mc2[1:firstn, 1]))
          comm_p = unique(unlist(mc2[,1]))
          biapattba_n = unique(unlist(mc2[1:firstn,1]))
          biapattba = unique(unlist(mc2[,1]))
        }
      }else{
        biapatt = append(biapatt, NA)
      }
    }
    comm_p_n = paste0(sort(comm_p_n), collapse = ', ')
    comm_p = paste0(sort(comm_p), collapse = ', ')
    biapattba_n = paste0(sort(biapattba_n), collapse = ', ')
    biapattba = paste0(sort(biapattba), collapse = ', ')
    new_frm[l,] = c(as.character(x2$Park_Name[1]), as.character(x2$FactorA[1]), as.character(x2$FactorB[1]),
                    as.character(x2$FactorC[1]), as.character(x2$FactorD[1]), x2$StationID[1], 
                    as.character(x2$VisitType[1]), x2$VisitId[1], x2$VisitDurMinutes[1], biapatt, uni_cod,
                    comm_p, comm_p_n, biapattba_n, biapattba)
  }
  return(new_frm)
}


#############################################################

### Code for Markov chain matrix for total time

mcmattime <- function(x, timevar, resvar, visit = 'VisitId',byvar = NULL){
  
  
  if (is.null(byvar)){
    
    dis_co = unlist(unique(x[resvar]))
    uni_visit = unlist(unique(x[visit]))
    len1 = length(uni_visit)
    len = length(dis_co)
    mcmatx = array(0, c(len, len), dimnames = list(dis_co, dis_co))
    for (j in 1:len1){
      x1 = x[x[visit] == uni_visit[j],]
      x1 = x1[order(x1[timevar], decreasing = FALSE),]
      n = dim(x1)
      if (n[1]>1){
        for (i in 2:n[1]){
          mcmatx[as.character(x1[resvar][i-1,]), as.character(x1[resvar][i,])] = mcmatx[as.character(x1[resvar][i-1,]), 
                                                                                        as.character(x1[resvar][i,])] + 
            as.numeric(difftime(x1[timevar][i,],x1[timevar][i-1,], units="hours"))
        }
      }
    }
  }else{
    dis_co = unlist(unique(x[resvar]))
    len = length(dis_co)
    bydis_co = unlist(unique(x[byvar]))
    len1 = length(bydis_co)
    mcmatx = array(0, c(len, len, len1), dimnames = list(dis_co, dis_co, bydis_co))
    for (j in 1:len1){
      x1 = x[x[byvar]== as.character(bydis_co[j]),]
      uni_visit = dis_co = unlist(unique(x1[visit]))
      len2 = length(uni_visit)
      for (k in 1:len2){
        x2 = x1[x1[visit]== uni_visit[k],]
        x2 = x2[order(x2[timevar], decreasing = FALSE),]
        n = dim(x2)
        if (n[1]>1){
          for (i in 2:n[1]){
            mcmatx[as.character(x2[resvar][i-1,]), 
                   as.character(x2[resvar][i,]), 
                   as.character(x2[byvar][i,])] = mcmatx[as.character(x2[resvar][i-1,]), 
                                                         as.character(x2[resvar][i,]), 
                                                         as.character(x2[byvar][i,])] + 
              as.numeric(difftime(x1[timevar][i,],x1[timevar][i-1,], units="hours"))
          }
        }
      }
    }
  }
  rm(dis_co, x, len, n, i)
  return(mcmatx)
  
}

tranc_tran_matrix <- function(x, trancval = 20){
  rsum = rowSums(x)
  y = x/rsum
  y[rsum > trancval,] = 0
  return(y)
}

# Transition Probability Matrix calculation

tran_prob_mat <- function(x){
  return(x/rowSums(x))
}

tran_prob_mat_all <- function(x){
  return(x/sum(x))
}

tran_prob_mat_lift <- function(x){
  return(x/(rowSums(x)*colSums(x)))
}
##########################################################


###################################################
## Code for creating Pattern Matrix

pattrn_mat <- function(x, varlist){
  require(reshape2)
  kk = data.frame(table(x[varlist]))
  name = names(kk)
  names(kk) = c('Var1', 'Var2', 'Freq')
  mcmatx = acast(kk, Var1 ~ Var2, value.var="Freq")
  # dis_co = unlist(unique(kk[,2]))
  # dis_visitID = unlist(unique(kk[,1]))
  # mcmatx = matrix(0, nrow = length(dis_visitID), ncol = length(dis_co))
  # rownames(mcmatx) = dis_visitID
  # colnames(mcmatx) = dis_co
  # 
  # for (i in 1:dim(kk)[1]){
  #   mcmatx[as.character(kk[i,1]), as.character(kk[i,2])] = kk[i,3]
  # }
  # rm(kk, dis_visitID, dis_co)
  return(mcmatx)
}

##################################################################

###################################################################
## Correlation of Pattern

cor_based_pattern <- function(x){
  return(paste0(names(x[order(x, decreasing = TRUE)[1:3]]), collapse = " -> "))
}

###################################################################
## subset best nmax probability

max_MC_erro <- function(x, nmax = 5){
  x1 = data.frame(as.table(prob_mat))
  x1 = x1[((x1$Freq != 0) & !is.nan(x1$Freq)),]
  d1 = dim(x1)[2]
  x1$xrank <- ave(-x1$Freq, x1[,c(d1-1,1:max(1,(d1-3)))], FUN=rank)
  #x1 = x1[order(x1$Var3,x1$Var1,-x1$Freq),]
  x1 = x1[x1$xrank %in% c(1:nmax),]
  
  return(x1)
}
################################################################
## Create a dataset with 5 most frequent pattern in order (Before, In, After)

pattern5seq <- function(data, len = c(1,2,2)){
  seq1 =seq(1:sum(len))
  bname = paste('BEr', seq1, sep = '')
  iname = paste('IEr', seq1, sep = '')
  aname = paste('AEr', seq1, sep = '')
  new_frm = data.frame(matrix(vector(), 0, sum(len)*3+13), stringsAsFactors = FALSE)
  names(new_frm) = c('Park_Name', 'FactorA', 'FactorB', 'FactorC', 'FactorD', 'StationID', 'VisitType', 
                     'VisitId', 'VisitDurMinutes' ,bname, iname, aname, 'Pattern','Time_Diff','Prev_TimeOn', 'Curr_TimeOn')
  l = 0
  aa = c('Before', 'In', 'After')
  unique_SID = unique(data$StationID)
  for (j in 1:length(unique_SID)){
    data11 = data[data$StationID == unique_SID[j],]
    data11 = data11[order(data11$Time_On),]
    unique_VisitId = unique(data11$VisitId)
    for (k in 1:length(unique_VisitId)){
      l = l + 1
      data1 = data11[data11$VisitId == unique_VisitId[k],]
      data1 = data1[order(data1$Time_On),]
      freq = data.frame(table(data1$Code_B_A, data1$EventWarningStop, data1$Code))
      freq = freq[order(freq[,1], freq[,2], -freq[,4]), ]
      bb = NULL
      
      for (i in aa){
        freq1 = freq[freq[,1]== i,]
        e1 = freq1[freq1[,2]=='Event' & freq1[,4]>0,][,3]
        w2 = freq1[freq1[,2]=='Warning' & freq1[,4]>0,][,3]
        s2 = freq1[freq1[,2]=='Stop' & freq1[,4]>0,][,3]
        if (i == 'Before'){
          if (length(e1) >= len[1]){
            e1 = paste('BE', e1[1:len[1]], sep = '')
          }else if (length(e1) > 0 & length(e1) < len[1]){
            e1 = paste('BE', e1[1:length(e1)], sep = '')
          }else{e1 = NA}
          if (length(w2) >= len[2]){
            w2 = paste('BW', w2[1:len[2]], sep = '')
          }else if (length(w2) >  0 & length(w2) < len[2]){ 
            w2 = paste('BW', w2[1:length(w2)], sep = '')
          }else{w2 = NA}
          if (length(s2) >= len[3]){
            s2 = paste('BS', s2[1:len[3]], sep = '')
          }else if (length(s2) > 0 & length(s2) < len[3]){ 
            s2 = paste('BS', s2[1:length(s2)], sep = '')
          }else{s2 = NA}
        }else if (i == 'In'){
          if (length(e1) >= len[1]){
            e1 = paste('IE', e1[1:len[1]], sep = '')
          }else if (length(e1) > 0 & length(e1) < len[1]){
            e1 = paste('IE', e1[1:length(e1)], sep = '')
          }else{e1 = NA}
          if (length(w2) >= len[2]){
            w2 = paste('IW', w2[1:len[2]], sep = '')
          }else if (length(w2) >  0 & length(w2) < len[2]){ 
            w2 = paste('IW', w2[1:length(w2)], sep = '')
          }else{w2 = NA}
          if (length(s2) >= len[3]){
            s2 = paste('IS', s2[1:len[3]], sep = '')
          }else if (length(s2) > 0 & length(s2) < len[3]){ 
            s2 = paste('IS', s2[1:length(s2)], sep = '')
          }else{s2 = NA}
        }else if (i == 'After'){
          if (length(e1) >= len[1]){
            e1 = paste('AE', e1[1:len[1]], sep = '')
          }else if (length(e1) > 0 & length(e1) < len[1]){
            e1 = paste('AE', e1[1:length(e1)], sep = '')
          }else{e1 = NA}
          if (length(w2) >= len[2]){
            w2 = paste('AW', w2[1:len[2]], sep = '')
          }else if (length(w2) >  0 & length(w2) < len[2]){ 
            w2 = paste('AW', w2[1:length(w2)], sep = '')
          }else{w2 = NA}
          if (length(s2) >= len[3]){
            s2 = paste('AS', s2[1:len[3]], sep = '')
          }else if (length(s2) > 0 & length(s2) < len[3]){ 
            s2 = paste('AS', s2[1:length(s2)], sep = '')
          }else{s2 = NA}
        }
        bb = append(bb, (c(c(e1[order(e1)], rep(NA,len[1]-length(e1))), c(w2[order(w2)], rep(NA,len[2]-length(w2))), 
                           c(s2[order(s2)], rep(NA,len[3]-length(s2))))))
      }
      if (k == 1){
        time_dif = 0
        time1 = as.Date(data1$Time_On_date[1], format = '%d/%m/%Y')
        time2 = as.Date(data1$Time_On_date[1], format = '%d/%m/%Y')
        time3 = c(format(time1, '%d/%m/%Y'), format(time2, '%d/%m/%Y'))
      }else {
        time2 = as.Date(data1$Time_On_date[1], format = '%d/%m/%Y')
        time_dif = as.numeric(time2-time1)
        time3 = c(format(time1, '%d/%m/%Y'), format(time2, '%d/%m/%Y'))
        time1 = time2
      }
      patt = paste(bb[!is.na(bb) & !(bb %in% "")], collapse = '.')
      new_frm[l,] = c(as.character(data1$Park_Name[1]), as.character(data1$FactorA[1]), as.character(data1$FactorB[1]),
                      as.character(data1$FactorC[1]), as.character(data1$FactorD[1]), data1$StationID[1], 
                      as.character(data1$VisitType[1]), data1$VisitId[1], data1$VisitDurMinutes[1], unname(bb), patt, 
                      time_dif, time3)
    }
  }
  return(new_frm)
}


#######################################################
# Cluster Analysis

clust <- function(data, nclust, byvar = NULL){
  
  if (is.null(byvar)){
    mm = data.frame(Code = numeric(0), Cluster_no = numeric(0));
    names(mm) = c('Code', 'Cluster_no')
    x = data
    for (i in 1:(nclust - 1)){
      aa = data.frame(table(x$VisitId, x$Code))
      a = acast(aa, Var1~Var2, value.var="Freq")
      a = scale(a)
      claster = kmeans((t(a)), 2, iter.max = 200, nstart = 100)
      b1 = table(claster$cluster)
      bb = data.frame(cbind(as.integer(names(claster$cluster)),claster$cluster))
      names(bb) = c('Code', 'Cluster_no')
      bb1 = bb[bb$Cluster_no == as.numeric(names(b1[b1 == min(b1)])),]
      bb1$Cluster_no = i
      mm = rbind(mm, bb1)
      x = x[!(x$Code %in% unique(mm$Code)),]
    }
    bb1 = bb[bb$Cluster_no == as.numeric(names(b1[b1 == max(b1)])),]
    bb1$Cluster_no = i+1
    mm = rbind(mm, bb1)
    return(mm)
  }else if (!is.null(byvar)){
    mm = data.frame(Code = numeric(0), Cluster_no = character(), by_var = character());
    names(mm) = c('Code', 'Cluster_no', byvar)
    uni_byvar = unlist(unique(data[byvar]))
    for (j in 1:length(uni_byvar)){
      x = data[data[byvar] == as.character(uni_byvar[j]),]
      for (i in 1:(nclust - 1)){
        aa = data.frame(table(x$VisitId, x$Code))
        a = acast(aa, Var1~Var2, value.var="Freq")
        a = scale(a)
        claster = kmeans((t(a)), 2, iter.max = 200, nstart = 100)
        b1 = table(claster$cluster)
        bb = data.frame(cbind(as.integer(names(claster$cluster)),claster$cluster))
        bb['by_var'] = uni_byvar[j]
        names(bb) = c('Code', 'Cluster_no', 'by_var')
        bb1 = bb[bb$Cluster_no == as.numeric(names(b1[b1 == min(b1)])),]
        bb1$Cluster_no = paste(uni_byvar[j], i)
        mm = rbind(mm, bb1)
        x = x[!(x$Code %in% unique(mm$Code)),]
      }
      bb1 = bb[bb$Cluster_no == as.numeric(names(b1[b1 == max(b1)])),]
      bb1$Cluster_no = paste(uni_byvar[j], i+1)
      mm = rbind(mm, bb1)
    }
    names(mm) = c('Code', 'Cluster_no', as.character(byvar))
    return(mm)
  }
}
