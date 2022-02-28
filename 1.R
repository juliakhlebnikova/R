library(ggplot2)
library(psych)
library(moments)

sunspot.year
write.csv2(sunspot.year, file="data.csv", row.names=T)

data <- sunspot.year

n1=sample(1:289, 50)
n2=sample(1:289, 20)

sub_data1=data[n1]
sub_data2=data[n2]

stat <- data.frame(matrix(NA,nrow=3, ncol=7))
colnames(stat) <- c("mean_value", "biased_variance", "unbiased_variance", 
                     "biased_deviation", "unbiased_deviation", "skewness", "kurtosis")


#Ex 1

k <- 1+round(3.322*log10(289))
d <- ceiling((max(data)-min(data))/k)


D=double(k)
D[1]=d-0.5
for(i in 2:k)
{ D[i]=D[i-1]+d }

Y=double(k)
Y[1]=0
for (i in 2:k) {Y[i]=D[i-1]}

data1 <- data.frame(interval=NA*c(1:k), middle_point=NA*c(1:k), frequency=NA*c(1:k), accumulated_frequency=NA*c(1:k))
sub1 <- data.frame(interval=NA*c(1:k), middle_point=NA*c(1:k), frequency=NA*c(1:k), accumulated_frequency=NA*c(1:k))
sub2 <- data.frame(interval=NA*c(1:k), middle_point=NA*c(1:k), frequency=NA*c(1:k), accumulated_frequency=NA*c(1:k))
data1$interval <- paste(Y,D,sep="-")
sub1$interval <- paste(Y,D,sep="-")
sub2$interval <- paste(Y,D,sep="-")

data1$middle_point[1] <- (Y[1] - 0.5 + D[1]) / 2
sub1$middle_point[1] <- (Y[1] - 0.5 + D[1]) / 2
sub2$middle_point[1] <- (Y[1] - 0.5 + D[1]) / 2
for (i in 2:k) {
  data1$middle_point[i] <- (Y[i] + D[i]) / 2
  sub1$middle_point[i] <- (Y[i] + D[i]) / 2
  sub2$middle_point[i] <- (Y[i] + D[i]) / 2
  # data1$middle_point[i] <- d/2 + (i-1)*d
}

for (i in 1:k) {
  data1$frequency[i] <- length(data[data >= ((i-1)*d) & data < (d*i)])
  sub1$frequency[i] <- length(sub_data1[sub_data1 >= ((i-1)*d) & sub_data1 < (d*i)])
  sub2$frequency[i] <- length(sub_data2[sub_data2 >= ((i-1)*d) & sub_data2 < (d*i)])
}

frequencies <- c(data1$frequency[1])
frequencies1 <- c(sub1$frequency[1])
frequencies2 <- c(sub2$frequency[1])
for ( i in 2:k ) {
  frequencies <- c(frequencies, frequencies[i-1]+data1$frequency[i])
  frequencies1 <- c(frequencies1, frequencies1[i-1]+sub1$frequency[i])
  frequencies2 <- c(frequencies2, frequencies2[i-1]+sub2$frequency[i])
  
}
data1["accumulated_frequency"] <- frequencies
sub1["accumulated_frequency"] <- frequencies1
sub2["accumulated_frequency"] <- frequencies2


# Ex 2

Fn <- ecdf(data)

plot(Fn, xlab="x", ylab="Fn(x)", main="Distribution Function", frame=FALSE)


#Ex 3

Fun <- function(arg) {
  n=double(1)
  for (i in 1:k)
  {
    if (arg < d*i & arg >= d*(i-1)) {n = data1$frequency / 289}
  }
  n=n/d
  return (n) 
} 

hist(data, breaks=k,freq=FALSE, main="Density Function", xlab="x", ylab="fn(x)", plot=TRUE, labels=FALSE)


# Ex 4

data = sort(data)
sub_data1 = sort(sub_data1)
sub_data2 = sort(sub_data2)

data_samples <- list(data, sub_data1, sub_data2)

for ( i in 1:3 ) {
  n <- length(data_samples[[i]])
  stat$mean_value[i] <- mean(data_samples[[i]])
  stat$biased_variance[i] <- sum((data_samples[[i]] - stat$mean_value[i]) ^ 2) / n
  stat$unbiased_variance[i] <- sum((data_samples[[i]] - stat$mean_value[i]) ^ 2) / (n - 1)
  stat$biased_deviation[i] <- sqrt(stat$biased_variance[i])
  stat$unbiased_deviation[i] <- sqrt(stat$unbiased_variance[i])
  stat$skewness[i] <- skewness(data_samples[[i]])
  stat$kurtosis[i] <- kurtosis(data_samples[[i]])
}


# Ex 5

mean_interval <- data.frame(matrix(NA, nrow=3, ncol=2))
colnames(mean_interval) <- c("left_bound", "right_bound")

for (i in 1:3) {
  n <- length(data_samples[[i]])
  mean_interval$left_bound[i] <- stat$mean_value[i] - qnorm (1-0.05/2)*stat$unbiased_deviation[i] / sqrt(n)
  mean_interval$right_bound[i] <- stat$mean_value[i] + qnorm (1-0.05/2)*stat$unbiased_deviation[i] / sqrt(n)
}

variance_interval <- data.frame(matrix(NA, nrow=3, ncol=2))
colnames(variance_interval) <- c("left_bound", "right_bound")

for (i in 1:3) {
  n <- length(data_samples[[i]])
  variance_interval$left_bound[i] <- (n-1) * stat$unbiased_variance[i] / qchisq(1-0.05/2, n-1)
  variance_interval$right_bound[i] <- (n-1) * stat$unbiased_variance[i] / qchisq(0.05/2, n-1)
}


# Ex 6

data_size <- data.frame(matrix(NA, nrow=3, ncol=2))
colnames(data_size) <- c("e1", "e2")

for (i in 1:3) {
  data_size$e1[i] <- (qnorm(1-0.05/2) * stat$unbiased_deviation[i]/0.001) ^ 2
  data_size$e2[i] <- (qnorm(1-0.05/2) * stat$unbiased_deviation[i]/0.0001) ^ 2
}

# Ex 7

# Pirson
suppose_data <- data1[c("middle_point", "frequency")]
l <- 1 / stat$mean_value[1]
suppose_data["p"] <- 0
suppose_data$p[1] <- 1 - exp(-l * (suppose_data$middle_point[1] + d / 2))
for (i in 2:k) {
  suppose_data$p[i] <- exp(-l * (suppose_data$middle_point[i] - d / 2)) - exp(-l * (suppose_data$middle_point[i] + d / 2))
}
suppose_data["pn"] <- sum(suppose_data$frequency) * suppose_data$p

result <- sum((suppose_data$freq - suppose_data$pn) ^ 2 / suppose_data$pn)

if (qchisq(1-0.05, k - 2) > result) {
  print("yes")
} else {
  print("no")
}


# Kolmogorov

suppose_data <- data1[c("middle_point", "frequency", "accumulated_frequency")]
n <- sum(suppose_data$frequency)

suppose_data["F"] <- 1 - exp(-l * suppose_data$middle_point)

suppose_data["Fn-F"] <- abs(suppose_data$accumulated_frequency / n - suppose_data$`F`)

lambda <- max(suppose_data$`Fn-F`) * sqrt(k)

# 1.36 for 1-0.05

if (lambda < 1.36) {
  print ("yes")
} else {
  print("no")
}


# Ex 8 
k1 <- sum(sub1$frequency)
k2 <- sum(sub2$frequency)
suppose_data <- sub2["middle_point"]
suppose_data[c("accumulated_frequency1", "accumulated_frequency2")] <- c(sub1$accumulated_frequency, sub2$accumulated_frequency)

suppose_data["diff"] <- abs(suppose_data$accumulated_frequency1 / k1 - suppose_data$accumulated_frequency2 / k2)
lambda <- sqrt(k1 * k2 / (k1 + k2)) * max(suppose_data$diff)

if (lambda < 1.36) {
  print("yes")
} else {
  print("no")
}


# Ex 9

suppose_data <- data1[c("middle_point", "frequency")]
n <- sum(suppose_data$frequency)

# exp val=45 lvl=0.01

z <- (stat$mean_value[1] - 45) / stat$unbiased_deviation[1] * sqrt(n)

if (z > qnorm(0.01) && z < qnorm(0.99)) {
  print("yes")
} else {
  print("no")
}


# Ex 10

# disp=1600 lvl=0.01

suppose_data <- data1[c("middle_point", "frequency")]
n <- sum(suppose_data$frequency)

z <- (n-1) * stat$unbiased_variance[1] / 1600

if (z < qchisq(0.99, n-1) && z > qchisq(0.01, n-1)) {
  print ("yes")
} else {
  print ("no")
}


# Ex 11


#мат ожид
n <- sum(data1$frequency)
n1 <- sum(sub1$frequency)
z <- (stat$mean_value[1] - stat$mean_value[2]) / sqrt(stat$unbiased_variance[1] / n + stat$unbiased_variance[2] / n1)
if (z > qt(0.05, n + n1 - 1) && z < qt(0.95, n + n1 - 1)) {
  print("yes")
} else {
  print("no")
}

#дисперсии, критерий фишера
z <- stat$unbiased_variance[2] / stat$unbiased_variance[1]
if (z > qf(0.05, n1-1, n-1) && z < qf(0.95, n1-1, n-1)) {
  print("yes")
} else {
  print("no")
}

###

