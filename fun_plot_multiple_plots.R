

# Draw multiple lines on the same plot ------------------------------------

# Set the seed for reproducibility
set.seed(123)

# Generate the data
x <- 1:25
y1 <- cumsum(rnorm(25))
y2 <- cumsum(rnorm(25))

# Create the plot
plot(x, y1, type = 'l', col = 'blue', ylim = c(min(y1, y2), max(y1, y2)), 
     xlab = 'X-axis', ylab = 'Y-axis', main = 'Overlaying Multiple Lines')
lines(x, y2, col = 'red')
legend('topleft', legend = c('Line 1', 'Line 2'), col = c('blue', 'red'), lty = 1)


# Side-by-Side Plots ------------------------------------------------------

# Create a side-by-side layout
par(mfrow = c(1, 2))

# Create the first plot
plot(x, y1, type = 'l', col = 'blue', 
     xlab = 'X-axis', ylab = 'Y-axis', main = 'Side-by-Side Plots (1)')

# Create the second plot
plot(x, y2, type = 'l', col = 'red',
     xlab = 'X-axis', ylab = 'Y-axis', main = 'Side-by-Side Plots (2)')

# Reset Par
par(mfrow = c(1, 1))


# Stacked Plots -----------------------------------------------------------

par(mfrow = c(2, 1), mar = c(2, 4, 4, 2))

# Create the first plot
plot(x, y1, type = 'l', col = 'blue', 
     xlab = 'X-axis', ylab = 'Y-axis', main = 'Stacked Plots')

# Create the second plot
plot(x, y2, type = 'l', col = 'red',
     xlab = 'X-axis', ylab = 'Y-axis', main = 'Stacked Plots (2)')

par(mfrow = c(1, 1))
