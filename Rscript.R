#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

#eins = args[3]	
cat("\nfirst argument:\t\n", args[1],"\n"  )
cat("second argument:\t\n", args[2],"\n")
cat("third argument:\t\n", args[3],"\n")


 
lst4R.PF_P <- read.delim(args[1]) 
lst4R.PF <- read.delim(args[2])
lst4R.PS <- read.delim(args[3])

#cat("lst4R.PF_P\n")
View(lst4R.PF_P)
View(lst4R.PF)
View(lst4R.PS)
cat("outcomes for the whole tree: \n\n")

sumPF_P_genes=sum(lst4R.PF_P$genes)
sumPF_genes=sum(lst4R.PF$genes)
sumPS_genes=sum(lst4R.PS$genes)

sumPF_P_gain=sum(lst4R.PF_P$gain)
sumPS_gain= sum(lst4R.PS$gain)
sumPF_gain=sum(lst4R.PF$gain)

sumPF_P_loss=sum(lst4R.PF_P$loss)
sumPS_loss=sum(lst4R.PS$loss)
sumPF_loss=sum(lst4R.PF$loss)

sumPF_P_gainloss=sum(lst4R.PF_P$gain)+sum(lst4R.PF_P$loss)
sumPF_gainloss=sum(lst4R.PF$gain)+sum(lst4R.PF$loss)
sumPS_gainloss=sum(lst4R.PS$gain)+sum(lst4R.PS$loss)

cat("total gain/loss events:\n\n")
cat("sumPF_P_genes:\t", sumPF_P_genes, "\n")
cat("sumPF_genes:\t", sumPF_genes, "\n")
cat("sumPS_genes:\t",sumPS_genes, "\n\n")

cat("sumPF_P_gain:\t",sumPF_P_gain,"\n")
cat("sumPF_gain:\t",sumPF_gain,"\n")
cat("sumPS_gain:\t",sumPS_gain,"\n")

cat("sumPF_P_loss:\t",sumPF_P_loss,"\n")
cat("sumPF_loss:\t",sumPF_loss,"\n")
cat("sumPS_loss:\t",sumPS_loss,"\n\n")

cat("sumPF_P_gainloss:\t",sumPF_P_gainloss,"\n\n")
cat("sumPF_gainloss:\t",sumPF_gainloss,"\n")
cat("sumPS_gainloss:\t",sumPS_gainloss,"\n")

#######plot total values:

gain=c(sumPS_gain,sumPF_gain,sumPF_P_gain)
names(gain)=c("PS","PF","PF_P")
loss=c(sumPS_loss,sumPF_loss,sumPF_P_loss)
gainloss=c(sumPS_gainloss,sumPF_gainloss,sumPF_P_gainloss)

total=cbind(gain,loss,gainloss)
names(total)=c("gain","loss","gain + loss")

path=paste(dirname(args[1]),"/totals_whole_tree.svg")

path2=gsub("[[:blank:]]", "", path)


#cat("path2",path2, "\n")

#pdf("~/epope_2_paper/epope_bash_run/results_for_bash_test/totals_whole_tree.pdf")
svg(path2)
mainhead="summation of total events\n-whole tree-\n"

par(mar=c(3.5,4.2,4.2,3), oma=c(0,0,0,0))
barplot(total, beside=T,  col = c("red","blue","#006600") , ylab = "summation of total events",
        main=mainhead,space=c(0.2,0.8),cex.lab=1.8, cex.axis=1.2,
        cex.main=1.8, cex.sub=1.8, cex.names=1.8,names.arg=c("gain", "loss", "gain + loss"))
legend(0.8,max(total),c("Sankoff", "Partition Function","Partition Function -P"),cex=1.3,
       fill=c("red","blue","#006600"))
abline(h=0)

dev.off()


#sum(lst4R.PS$Nr_Fams)

#sum(lst4R.PF$Nr_Fams)
#sum(lst4R.PF_P$Nr_Fams)

#sum(lst4R.PS$gainFam)
#sum(lst4R.PF$gainFam)
#sum(lst4R.PF_P$gainFam)

#sum(lst4R.PS$lossFam)
#sum(lst4R.PF$lossFam)
#sum(lst4R.PF_P$lossFam)

calculate = function(x)
{ if (class(x) != "data.frame")
  stop("this should be a dataframe")
  #else print("allright")
  len = length(x$lab)
  len=c(1:len)
  x$species=as.character(x$species)
  #suche den vater, wenn vater nicht gefunden, berechne nächsten
  for (i in len)
  {
    tempparlab=x$parents[i]
    found=0
        #cat ("vater von ", x$species[i], "ist", tempparlab, "\n" )
    # diese parentslab müssen jetzt in lab gefunden werden
    for (j in len)
    { if(tempparlab == x$lab[j]) 
    { found=1
      #cat ("Vater von:", x$species[i], "ist ", x$species[j], "\n")
      #hier kann man alles ausrechnen, i hält gain loss zeug und j den vater davon:
      x$gain[i]   =(100/x$genes[i])*x$gain[i]
     #cat ("gain: ", x$gain[i], "species:",x$species[i],"\n")
      if(is.nan(x$gain[i]))
      {
        x$gain[i]=0
      }
      x$loss[i]   =(100/x$genes[j])*x$loss[i]
      if(is.nan(x$loss[i]))
      {
        x$loss[i]=0
      }
      x$gainFam[i]=(100/x$Nr_Fams[i])*x$gainFam[i]
      if(is.nan(x$gainFam[i]))
      {
        x$gainFam[i]=0
      }
      x$lossFam[i]=(100/x$Nr_Fams[j])*x$lossFam[i]
      if(is.nan(x$lossFam[i]))
      {
        x$lossFam[i]=0
      }
    }
    }
    if (found == 0)
      {
        #cat("no parent for ",x$species[i] ,"\n")
        if(x$gain[i]>100)
          x$gain[i]=100
        
        
      }
      
  }
  return (x)
}

lst4R.PS=calculate(lst4R.PS)
lst4R.PF_P=calculate(lst4R.PF_P)
lst4R.PF=calculate(lst4R.PF)
#View(lst4R.PF)

cat("sum of relative gain/loss events:\n\n")

cat("sumPS_gain:\t",sum(lst4R.PS$gain ),"\n")
cat("sumPF_gain:\t",sum(lst4R.PF$gain ),"\n")
cat("sumPF_P_gain:\t",sum(lst4R.PF_P$gain),"\n\n")

cat("sumPS_loss:\t",sum(lst4R.PS$loss),"\n")
cat("sumPF_loss:\t",sum(lst4R.PF$loss),"\n")
cat("sumPF_P_loss:\t",sum(lst4R.PF_P$loss),"\n\n")

cat("sumPS_gainloss:\t",sum(lst4R.PS$loss)+sum(lst4R.PS$gain ),"\n")
cat("sumPF_gainloss:\t",sum(lst4R.PF$gain)+sum(lst4R.PF$loss),"\n")
cat("sumPF_P_gainloss:\t",sum(lst4R.PF_P$gain)+sum(lst4R.PF_P$loss),"\n\n")

#~ #######plot relative  values:

path=paste(dirname(args[1]),"/relatives_whole_tree.svg")

path2=gsub("[[:blank:]]", "", path)


#pdf("~/epope_2_paper/epope_bash_run/results_for_bash_test/relatives_whole_tree.pdf")
svg(path2)
rel_gain=c(sum(lst4R.PS$gain ),sum(lst4R.PF$gain ),sum(lst4R.PF_P$gain))
names(rel_gain)=c("PS","PF","PF_P")

rel_loss=c(sum(lst4R.PS$loss),sum(lst4R.PF$loss),sum(lst4R.PF_P$loss))
rel_gainloss=c(sum(lst4R.PS$loss)+sum(lst4R.PS$gain ),sum(lst4R.PF$gain)+sum(lst4R.PF$loss),sum(lst4R.PF_P$gain)+sum(lst4R.PF_P$loss))
relative=cbind(rel_gain,rel_loss,rel_gainloss)
relative
mainhead="summation of relative events\n-whole tree-\n"

par(mar=c(3.5,4.2,4.2,3), oma=c(0,0,0,0))
barplot(relative, beside=T,col = c("red","blue","#006600") , ylab = "summation of relative events",
        main=mainhead,space=c(0.2,0.8),cex.lab=1.8, cex.axis=1.2, cex.main=1.8, cex.sub=1.8, cex.names=1.8,
        names.arg = c("gain","loss","gain + loss"))
legend(0.8,max(relative),c("Sankoff", "Partition Function" , 
                           "Partition Function -P"),fill=c("red","blue","#006600"),cex=1.3)
abline(h=0)
dev.off()

cat("absolute events: \n\n")
cat("PF_loss > than PF_P_loss:",sum(lst4R.PF$loss>lst4R.PF_P$loss),"\n")
cat("PF_loss < than PF_P_loss:", sum(lst4R.PF$loss<lst4R.PF_P$loss),"\n")
cat("PF_loss == than PF_P_loss:", sum(lst4R.PF$loss==lst4R.PF_P$loss),"\n\n")

cat("PF_gain > than PF_P_gain:", sum(lst4R.PF$gain>lst4R.PF_P$gain),"\n")
cat("PF_gain < than PF_P_gain:", sum(lst4R.PF$gain<lst4R.PF_P$gain),"\n")
cat("PF_gain == than PF_P_gain:", sum(lst4R.PF$gain==lst4R.PF_P$gain), "\n")

PS_PF_gain = data.frame(lst4R.PS$species,lst4R.PS$gain,lst4R.PF$gain,lst4R.PF_P$gain)
names(PS_PF_gain) = c("species","gain_PS","gain_PF","gain_PF_P")
PS_PF_gain 
# sum(PS_PF_gain$gain_PF)
# sum(PS_PF_gain$gain_PF_P)

PS_PF_loss = data.frame(lst4R.PF$species,lst4R.PS$loss,lst4R.PF$loss,lst4R.PF_P$loss)
names(PS_PF_loss) = c("species","loss_PS","loss_PF","loss_PF_P")

#View(PS_PF_loss)
PS_PF_gainFam = data.frame(lst4R.PS$species,lst4R.PS$gainFam,lst4R.PF$gainFam,lst4R.PF_P$gainFam)
names(PS_PF_gainFam) = c("species","gainFam_PS","gainFam_PF","gain_Fam_PF_P")

PS_PF_lossFam = data.frame(lst4R.PF$species,lst4R.PS$lossFam,lst4R.PF$lossFam,lst4R.PF_P$lossFam)
names(PS_PF_lossFam) = c("species","lossFam_PS","lossFam_PF","lossFam_PF_P")

#loss einfach negativ:
PS_PF_loss[,c(2,3,4)]=PS_PF_loss[,c(2,3,4)]* -1
PS_PF_loss
sum(PS_PF_loss$loss_PF)
sum(PS_PF_loss$loss_PF_P)
PS_PF_lossFam[,c(2,3,4)]=PS_PF_lossFam[,c(2,3,4)]* -1
PS_PF_lossFam

####################plot loss and gain:

plot_gain_PS=PS_PF_gain$gain_PS    #herausziehen der spalte gain von PS
plot_gain_PF=PS_PF_gain$gain_PF #herausziehen der spalte gain von PF
plot_gain_PF_P=PS_PF_gain$gain_PF_P
names(plot_gain_PF)=PS_PF_gain$species #werte mit speziesnamen benennen
names(plot_gain_PS)=PS_PF_gain$species #werte mit speziesnamen benennen
names(plot_gain_PF_P)=PS_PF_gain$species
grafik_gain_3=rbind(plot_gain_PF_P, plot_gain_PF,plot_gain_PS) #3 werte unter einem species namen jeweils zusammenführen
grafik_gain_3

plot_loss_PS=PS_PF_loss$loss_PS
plot_loss_PF=PS_PF_loss$loss_PF
plot_loss_PF_P=PS_PF_loss$loss_PF_P
names(plot_loss_PS)=PS_PF_loss$species 
names(plot_loss_PF)=PS_PF_loss$species  
names(plot_loss_PF_P)=PS_PF_loss$species
grafik_loss_3=rbind(plot_loss_PF_P,plot_loss_PF,plot_loss_PS)
grafik_loss_3

class(grafik_gain_3)
grafik_gain_3

#install.packages("extrafont")

#library(extrafont)
#font_import()

path=paste(dirname(args[1]),"/gain_loss_whole_tree.svg")

path2=gsub("[[:blank:]]", "", path)


  #pdf("~/epope_2_paper/epope_bash_run/results_for_bash_test/gain_loss_P_whole_tree.pdf",height=7,width = 31)
#für kleinen baum:
#svg("~/epope_2_paper/epope_bash_run/results_for_bash_test/gain_loss_P_whole_tree.svg",height=15 ,width =8.2  )
#für ganzen baum: 
svg(path2,height=32 ,width =5  )


#par(mar=c(5,4,4,4), oma=c(5,0,0,0))

#für ganzen baum: 
#par(mar=c(5,4,4,4), oma=c(5,4,4,4),cex=0.6)
#für kleinen baum
#par(mar=c(5,4,4,4), oma=c(5,4,4,4),cex=1.4,font=1)#font für legend
#für tunicaten baum
par(mar=c(5,4,4,4), oma=c(5,4,4,4),cex=0.6,font=1)#font für legend

barplot( grafik_gain_3,horiz=T, beside=T, xlim=c(-100,100),border=F, las=1,col = c("#006600","blue","red"))


mainhead2="relative gain and loss" #for parsimony, partition function\n and partition function with -P option\n -whole tree-"

barplot(grafik_loss_3, horiz=T, beside=T, las=1,add=T,border=F,xaxt="n",yaxt="n",
        col = c("#006600","blue","red") , 
        #legend=c("parsimony score", "partition function", "partition function -P"),
        xlab = "loss and gain in %",
        main=mainhead2)
legend("bottomright",legend=c("Sankoff", "Partition Function", "Partition Function -P"),col = c("red","blue","#006600")
        ,fill = c("red","blue","#006600"),cex = 1.1)

#für kleinen baum: cex=0.65
#für ganzen baum: cex=1.1
abline(v=0)
dev.off()





 #~ ############################plot loss and gain FAM

#~ # 
#~ # plot_gainFam_PS=PS_PF_gainFam$gainFam_PS    #herausziehen der spalte gain von PS
#~ # plot_gainFam_PF=PS_PF_gainFam$gainFam_PF #herausziehen der spalte gain von PF
#~ # names(plot_gainFam_PF)=PS_PF_gainFam$species #werte mit speziesnamen benennen
#~ # names(plot_gainFam_PS)=PS_PF_gainFam$species #werte mit speziesnamen benennen
#~ # grafik_gainFam_both=rbind(plot_gainFam_PS, plot_gainFam_PF) #beide werte unter einem species namen jeweils zusammenführen
#~ # grafik_gainFam_both
#~ # 
#~ # plot_lossFam_PS=PS_PF_lossFam$lossFam_PS
#~ # plot_lossFam_PF=PS_PF_lossFam$lossFam_PF
#~ # names(plot_lossFam_PS)=PS_PF_lossFam$species 
#~ # names(plot_lossFam_PF)=PS_PF_lossFam$species  
#~ # grafik_lossFam_both=rbind(plot_lossFam_PS, plot_lossFam_PF)
#~ # grafik_lossFam_both
#~ # 
#~ # if(wholetree){pdf("~/Dropbox/Masterarbeit/daten/NEUNEUNEU/both/gain_lossFam.printorder.pdf", height = 7, width = 31)}
#~ # if(tun){pdf("~/Dropbox/Masterarbeit/daten/NEUNEUNEU/both/gain_lossFam.tun.pdf",height=7,width=15 )}
#~ # if(gnat){pdf("~/Dropbox/Masterarbeit/daten/NEUNEUNEU/both/gain_lossFam.gnath.pdf", height=7, width=15)}
#~ # if(nem){pdf("~/Dropbox/Masterarbeit/daten/NEUNEUNEU/both/gain_lossFam.nem.pdf")}
#~ # 
#~ # par(mar=c(5,4,4,2), oma=c(6,0,0,0))
#~ # 
#~ # barplot(grafik_gainFam_both, beside=T, ylim=c(-100,100),border=F,
#~ #         las=3, 
#~ #         col = c("red","blue"))
#~ # 
#~ # if(wholetree)
#~ # {mainhead="gainFam and lossFam for parsimony score \n and partition function - whole tree "}
#~ # if(tun)
#~ # {mainhead="gainFam and lossFam for parsimony score \n and partition function - tunicata subtree "}
#~ # if(nem)
#~ # {mainhead="gainFam and lossFam for parsimony score \n and partition function - nematoda subtree "}
#~ # if(gnat)
#~ # {mainhead="gainFam and lossFam for parsimony score \n and partition function - gnathostomata subtree"}
#~ # 
#~ # barplot(grafik_lossFam_both, beside=T, las=3,add=T,border=F,legend= c("parsimony score", "partition function"),
#~ #         col = c("red","blue"),ylab = "gain and loss in %",main=mainhead)
#~ # 
#~ # abline(h=0)
#~ # dev.off()
#~ # 
#~ # 


