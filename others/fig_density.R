rm(list = ls())


# set path as current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

ngrid = 1000
xgrid = seq(0.0001, 0.9999, length = ngrid)

mygroup = as.factor(rep(c("cobin","micobin","beta"), each = ngrid))
# change order of levels
mygroup = factor(mygroup, levels = c("micobin","cobin","beta"))
library(ggplot2)
library(cobin)
library(betareg)
# cobin with lambda = 1 (conti Bernoulli), where mean = 0.3 and mean = 0.8
lambda1 = 1
mu1 = 0.3;
#eta1 = g(mu1)
eta1 = bftprimeinv(mu1)
var1 = cobin:::vcobin(eta1, lambda1)
# https://www.johndcook.com/blog/2021/04/07/beta-given-mean-variance/; phi is a+b
phi1 = mu1*(mu1*(1-mu1)/var1 -1) + (1-mu1)*(mu1*(1-mu1)/var1 -1)


plot(xgrid, dcobin(xgrid, eta1, lambda1), type = "l", col = 1, ylim = c(0,5))
lines(xgrid, dmicobin(xgrid, eta1, psi = 1/lambda1), col  = 2, lty = 2)
lines(xgrid, dbetar(xgrid, mu = mu1, phi1), col  = 3, lty = 3)

# plot(xgrid, dcobin(xgrid, eta1, lambda1, log = T), type = "l", col = 1, ylim = c(-20,2))
# lines(xgrid, dmicobin(xgrid, eta1, psi = 1/lambda1, log = T), col  = 2, lty = 2)
# lines(xgrid, dbetar(xgrid, mu = mu1, phi1, log= T), col  = 3, lty = 3)



a1 = mu1*phi1
b1 = (1-mu1)*phi1

round(eta1,2)


round(a1, 2)
round(b1, 2)
eta1
lambda1

mu1
round(var1,3)

df1 = data.frame(x = rep(xgrid,3),
                group = mygroup,
                #linetype = rep(c("soild","longdash","twodash"), each = 200),
                y = c(dcobin(xgrid, eta1, lambda1),
                      dmicobin(xgrid, eta1, psi = lambda1),
                      dbetar(xgrid, mu = mu1, phi1)))
mainplot1 = ggplot(df1, aes(x=x, y = y)) +
  geom_line(aes(col = group, linetype = group)) +
  xlab("y") + ylab("density") + ylim(0,4)+
  facet_wrap(~"mean = 0.3, var = 0.06") +
  scale_color_manual(labels =  c("cobin"="cobin(-2.67,1)", "micobin"="micobin(-2.67, 1)",  "beta"="beta(1.74,0.74)"),
                     values = c("cobin" = "#fc8d59", "micobin" = "black","beta"="#619CFF"),
                     breaks = c("beta","cobin","micobin"))+
  scale_linetype_manual(labels =  c("cobin"="cobin(-2.67,1)", "micobin"="micobin(-2.67, 1)",  "beta"="beta(1.74,0.74)"),
                        values = c("cobin" = 5, "micobin" = 1,"beta"=4),
                        breaks = c("beta","cobin","micobin")) +
  theme_bw() +
  theme(legend.title=element_blank(), legend.key.width = unit(1.5, "line")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7,0.8)) + theme(legend.text=element_text(size=7))+ theme(axis.title.x=element_blank())
  #theme(axis.title.x=element_blank())

mainplot1

lambda2 = 3
mu2 = 0.5;


#eta2 = g(mu2)
eta2 = bftprimeinv(mu2)
var2 = cobin:::vcobin(eta2, lambda2)
mu2
var2
# https://www.johndcook.com/blog/2021/04/07/beta-given-mean-variance/; phi is a+b
phi2 = mu2*(mu2*(1-mu2)/var2 -1) + (1-mu2)*(mu2*(1-mu2)/var2 -1)


a2 = mu2*phi2
b2 = (1-mu2)*phi2

round(eta2,2)


round(a2, 2)
round(b2, 2)


plot(xgrid, dcobin(xgrid, eta2, lambda2), type = "l", col = 1, ylim = c(0,5))
lines(xgrid, dmicobin(xgrid, eta2, psi = 1/lambda2), col  = 2, lty = 2)
lines(xgrid, dbetar(xgrid, mu = mu2, phi2), col  = 3, lty = 3)


df2 = data.frame(x = rep(xgrid,3),
                 group = mygroup,
                 #linetype = rep(c("soild","longdash","twodash"), each = 200),
                 y = c(dcobin(xgrid, eta2, lambda2),
                       dmicobin(xgrid, eta2, psi = 1/lambda2),
                       dbetar(xgrid, mu = mu2, phi2)))
mainplot2 = ggplot(df2, aes(x=x, y = y)) +
  geom_line(aes(col = group, linetype = group)) +
  xlab("y") + ylab("density") + ylim(0,4)+
  facet_wrap(~"mean = 0.5, var = 0.028") +
  scale_color_manual(labels =  c("cobin"="cobin(0,1/3)", "micobin"="micobin(0, 1/3)",  "beta"="beta(4,4)"),
                     values = c("cobin" = "#fc8d59", "micobin" = "black","beta"="#619CFF"),
                     breaks = c("beta","cobin","micobin"))+
  scale_linetype_manual(labels =  c("cobin"="cobin(0,1/3)", "micobin"="micobin(0, 1/3)",  "beta"="beta(4,4)"),
                        values = c("cobin" = 5, "micobin" = 1,"beta"=4),
                        breaks = c("beta","cobin","micobin")) +
  theme_bw() +
  theme(legend.title=element_blank(), legend.key.width = unit(1.5, "line")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.5,0.85)) + theme(legend.text=element_text(size=7)) + guides(colour = guide_legend(nrow = 2)) + theme(axis.title.x=element_blank())
  #theme(axis.title.x=element_blank())
mainplot2
#
# df3 = data.frame(x = rep(xgrid,3),
#                  group = rep(c("cobin","micobin","beta"), each = ngrid),
#                  #linetype = rep(c("soild","longdash","twodash"), each = 200),
#                  y = c(dcobin(xgrid, eta2, lambda2 ,log = T),
#                        dmicobin(xgrid, eta2, psi = 1/lambda2, log = T),
#                        dbetar(xgrid, mu = mu2, phi2, log = T)))
# mainplot3 = ggplot(df3, aes(x=x, y = y)) +
#   geom_line(aes(col = group, linetype = group)) +
#   xlab("") + ylab("log density") + theme_light() +
#   facet_wrap(~"mean = 0.5, var = 0.028") +
#   scale_color_manual(labels = c("cobin"="continuous binomial (cobin)", "micobin"="mixture of continuous binomial (micobin)",  "beta"="beta"),
#                      values = c("cobin" = "#fc8d59", "micobin" = "black","beta"="#99d594"))+
#   scale_linetype_manual(labels = c("cobin"="continuous binomial (cobin)", "micobin"="mixture of continuous binomial (micobin)",  "beta"="beta"),
#                         values = c("cobin" = 5, "micobin" = 4,"beta"=1)) +
#   theme(legend.title=element_blank(), legend.key.width = unit(1.5, "line")) +
#   theme(axis.title.x=element_blank())
# mainplot3

#
# lambda3 = 4
# mu3 = 0.7;

lambda3 = 2
mu3 = 0.8;


eta3 = bftprimeinv(mu3)
#eta3 = g(mu3)
var3 = cobin:::vcobin(eta3, lambda3)
mu3
var3
# https://www.johndcook.com/blog/2021/04/07/beta-given-mean-variance/; phi is a+b
phi3 = mu3*(mu3*(1-mu3)/var3 -1) + (1-mu3)*(mu3*(1-mu3)/var3 -1)


plot(xgrid, dcobin(xgrid, eta3, lambda3), type = "l", col = 1, ylim = c(0,5))
lines(xgrid, dmicobin(xgrid, eta3, psi = 1/lambda3), col  = 2, lty = 2)
lines(xgrid, dbetar(xgrid, mu = mu3, phi3), col  = 3, lty = 3)

a3 = mu3*phi3
b3 = (1-mu3)*phi3

round(eta3,2)


round(a3, 2)
round(b3, 2)

df3 = data.frame(x = rep(xgrid,3),
                 group = mygroup,
                 #linetype = rep(c("soild","longdash","twodash"), each = 200),
                 y = c(dcobin(xgrid, eta3, lambda3),
                       dmicobin(xgrid, eta3, psi = 1/lambda3),
                       dbetar(xgrid, mu = mu3, phi3)))
mainplot3 = ggplot(df3, aes(x=x, y = y)) +
  geom_line(aes(col = group, linetype = group)) +
  xlab("y") + ylab("density") + ylim(0,4)+
  facet_wrap(~"mean = 0.8, var = 0.018") +
  scale_color_manual(labels =  c("cobin"="cobin(4.80,1/2)", "micobin"="micobin(4.80, 1/2)",  "beta"="beta(6.51,1.63)"),
                     values = c("cobin" = "#fc8d59", "micobin" = "black","beta"="#619CFF"),
                     breaks = c("beta","cobin","micobin"))+
  scale_linetype_manual(labels =  c("cobin"="cobin(4.80,1/2)", "micobin"="micobin(4.80, 1/2)",  "beta"="beta(6.51,1.63)"),
                        values = c("cobin" = 5, "micobin" = 1,"beta"=4),
                        breaks = c("beta","cobin","micobin")) +
  theme_bw() +
  theme(legend.title=element_blank(), legend.key.width = unit(1.5, "line")) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.3,0.8)) + theme(legend.text=element_text(size=7)) + theme(axis.title.x=element_blank())
mainplot3




ggpubr::ggarrange(mainplot1, mainplot2,  mainplot3, ncol = 3)



#ggsave("fig1_new.pdf", width =  9, height = 2.8)



