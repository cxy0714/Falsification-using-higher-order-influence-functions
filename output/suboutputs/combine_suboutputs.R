load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,1281.RData')
sim_out1 <- sim_out
load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,1282.RData')
sim_out2 <- sim_out
load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,1283.RData')
sim_out3 <- sim_out
load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,1284.RData')
sim_out4 <- sim_out
load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,1285.RData')
sim_out5 <- sim_out
load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,1286.RData')
sim_out6 <- sim_out
load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,1287.RData')
sim_out7 <- sim_out
load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,1288.RData')
sim_out8 <- sim_out
load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,1289.RData')
sim_out9 <- sim_out
load('./output/suboutputs/continuous_sim_5000_20_2,8,32,64,12810.RData')
sim_out10 <- sim_out

sim_out <- cbind(sim_out1,sim_out2,sim_out3,sim_out4,sim_out5,
                 sim_out6,sim_out7,sim_out8,sim_out9,sim_out10) 

save(sim_out, file="./output/continuous_sim_5000_200_2,8,32,64,128.RData")


load('./output/suboutputs/continuous_sim_choose_stopping_k_5000_50_2,8,32,64,1281.RData')
sim_out1 <- sim_out
load('./output/suboutputs/continuous_sim_choose_stopping_k_5000_50_2,8,32,64,1282.RData')
sim_out2 <- sim_out
load('./output/suboutputs/continuous_sim_choose_stopping_k_5000_50_2,8,32,64,1283.RData')
sim_out3 <- sim_out
load('./output/suboutputs/continuous_sim_choose_stopping_k_5000_50_2,8,32,64,1284.RData')
sim_out4 <- sim_out

sim_out <- cbind(sim_out1,sim_out2,sim_out3,sim_out4) 

save(sim_out, file="./output/continuous_sim_choose_stopping_k_5000_200_2,8,32,64,128.RData")