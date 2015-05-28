### Batch adjust and analyse random data
#
# Compare F statistics obtained from ANOVA against:
# - reference: distribution assumed by common ANOVA
# - design corrected: F scaled to correct for effective sample size
# - design corrected (shrinked):F scaled assuming shrinked batch adjustment
# - empirically corrected: F scaled based on empirical factor (by variance)
# - theoretical value: approximate F distribution after batch adjustment


# Set work directory (needs to by modified!)
setwd('D:/Store/Git/batch-adjust-warning-reports/theory');

# Include libraries and functinos
library(limma);
library(MASS);
library(sva); # Not in use unless ComBat is used
source('../commonscripts/theoryfunctions.r');

# Settings
options(digits=5);

### Specify random data (N=design, dt=data)
#N = UnbalancedBatches(5,50,10,full=0,multi=1);
#N = list(c(10,10,0,0,0),c(0,10,10,0,0),c(0,0,10,10,0),c(0,0,0,10,10));
N = list(c(10,10,0,0,0),c(0,10,10,0,0),c(0,0,10,10,0),c(0,0,0,10,10),c(10,0,0,0,10));
dt = RandomData(N,1000,sd.batch=.5,sd.error=1);

### Batch adjustment
batch = dt$pheno$Batch;
group = dt$pheno$Group;
mod = model.matrix(~1+Group,data=dt$pheno); # Model (group effect)
mod0 = model.matrix(~1,data=dt$pheno); # Null-model (intercept)
mod.batch = model.matrix(~Batch,data=dt$pheno); # Batch effect model
dt$edata.anova = removeBatchEffect(dt$edata, batch=batch, design=mod);
#dt$ComBat(mod); # Run ComBat and store in dt.edata.combat

### ANOVA
dt$edata=dt$edata.anova;
#dt$edata=dt$edata.combat;
mod1 = mod[,2:ncol(mod)];
test = Ftest(dt$edata,mod,mod0);
B.A=t(mod1)%*%perpendicular(mod0)%*%mod1;
B.AC=t(mod1)%*%perpendicular(mod.batch)%*%mod1;
M=B.A%*%solve(B.AC);
Meig=eigen(M,only.values=TRUE)$values;
M1.est=tr(M); M2.est=tr(M%*%M);
df1.est=M1.est**2/M2.est; var.est=M2.est/M1.est;

### QQ-plots
test$N.adjust.emp=mean(test$MS0)/mean(test$MS1);
test$df0.adjust=dt$N.samples-dt$N.groups-dt$N.batches+1;
test$df1.adjust=df1.est;
ratio.est=M1.est/(dt$N.groups-1);
#test$ratio.est=var.est*(1+(dt$N.batches-1)/test$df0.adjust);
#Fquant=qf(ppoints(dt$N.genes),test$df1,test$df0);# Standard F-distribution
Fquant=qf(ppoints(dt$N.genes),test$df1.adjust,test$df0.adjust);# Adjusted F-distribution
qqplot(sqrt(Fquant),sqrt(test$F),col='blue');
Line.reset();
Line.add(0,1,col='gray',lwd=2,label='reference (diagonal)');
Line.add(0,1/sqrt(dt$N.adjust),col='red',lwd=2,label='design corrected (shrinked)');
Line.add(0,1/sqrt(dt$N.adjust0),col='orange',lwd=1,label='design corrected');
Line.add(0,1/sqrt(test$N.adjust.emp),col='skyblue',lwd=2,label='empirically corrected');
Line.add(0,sqrt(ratio.est),col='green',lwd=2,label='theoretical value');
legend('topleft',cex=.8,legend=Line.leg,lwd=2,col=Line.col);

M1.emp=mean(test$SS1);M2.emp=var(test$SS1)/2;
df1.emp=M1.emp**2/M2.emp; var.emp=M2.emp/M1.emp; ratio.emp=1/test$N.adjust.emp;
print.vars('ESTIMATES FROM THEORY:',M1.est,M2.est,df1.est,var.est,ratio.est);
print.vars('EMPIRICAL ESTIMATES:',M1.emp,M2.emp,df1.emp,var.emp,ratio.emp);
print.vars('EIGENVALUES OF M:',Meig);

