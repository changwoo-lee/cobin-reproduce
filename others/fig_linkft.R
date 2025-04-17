rm(list = ls())


# set path as current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(cobin)
library(betareg)

xgrid = seq(-20, 20, length = 500)

df1 = data.frame(x = rep(xgrid,3),
                group = rep(c("cobin","cauchit","logit"), each = 500),
                #linetype = rep(c("soild","longdash","twodash"), each = 200),
                y = c(bftprime(xgrid),
                      atan(xgrid*pi/12)/pi + 0.5,
                      1/(1+exp(-xgrid/3))))
mainplot1 = ggplot(df1, aes(x=x, y = y)) +
  geom_line(aes(col = group, linetype = group)) +
  xlab(expression(eta)) +
  ylab("") +
  theme_bw() +
 # theme(legend.justification = c(0, 1),
#        legend.position = c(0.05,0.95))+
  theme(legend.title=element_blank(), legend.key.width = unit(1.5, "line")) +
  facet_wrap(~"Comparison of inverse of link functions") +
  scale_color_manual(
    values = c("cobin" = "black", "cauchit" = "red", "logit" = "#619CFF"),
    breaks = c("cobin", "cauchit", "logit"),
    labels = c(
      "logit" = expression(1/(1+exp(-x/3))),
      "cobin" = expression(paste("B'(",x,") = ",exp(x)/(exp(x)-1)-1/x)),
      "cauchit" = expression(arctan((pi/12)*x)/pi + 0.5)
    )
  ) +
  scale_linetype_manual(
    values = c("cobin" = 1, "cauchit" = 3, "logit" = 4),
    breaks = c("cobin", "cauchit", "logit"),
    labels = c(
      "logit" = expression(1/(1+exp(-x/3))),
      "cobin" = expression(paste("B'(",x,") = ",exp(x)/(exp(x)-1)-1/x)),
      "cauchit" = expression(arctan((pi/12)*x)/pi + 0.5)
    )
  )+
  xlab(expression(x))+
  theme(legend.position = "top")
  #
mainplot1

xgrid2 = seq(0.001, 0.999, length = 500)

df2 = data.frame(x = rep(xgrid2,2),
                 group = rep(c("cobin","logit"), each = 500),
                 #linetype = rep(c("soild","longdash","twodash"), each = 200),
                 y = c(Vft(xgrid2),
                       xgrid2*(1-xgrid2)/3))

mainplot2 = ggplot(df2, aes(x=x, y = y)) +
  geom_line(aes(col = group, linetype = group)) +
  xlab(expression(mu)) +
  ylab("") +  theme_bw() +
facet_wrap(~"Comparison of variance functions") +
  theme(legend.title=element_blank(), legend.key.width = unit(1.5, "line")) +
  scale_color_manual(labels = c("cobin" = expression(V(mu)),
                                "logit" = expression(mu(1-mu)/3)),
                     values = c("cobin" = "black", "logit"="#619CFF"),
                     breaks = c("cobin",  "logit")) + # desired order)+
  scale_linetype_manual(labels = c("cobin" = expression(V(mu)),
                                   "logit" = expression(mu(1-mu)/3)),
                        values = c("cobin" = 1, "logit"=4),
                        breaks = c("cobin", "logit")) +

  xlab(expression(mu))+
  theme(legend.position = "top")

mainplot2




ggpubr::ggarrange(mainplot1, mainplot2, ncol = 2, widths = c(5, 3))

#ggsave("fig2_new.pdf", width =  10, height = 3)



