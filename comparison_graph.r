library(ggplot2)

summary_file <- read.csv('/Users/vikas/Downloads/bootcamp/metagenomics/bmebootcamp-metagenomics/read_comparison.csv', stringsAsFactors = FALSE)
summary_file$read_counts<-as.numeric(summary_file$read_counts)

plot<-ggplot(summary_file, aes(x = organism, y = read_counts, group = method, fill = method)) + 
  geom_histogram(stat = 'identity', position = 'dodge') + 
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(axis.text.x = element_text(angle = 45 , hjust = 1), axis.title.x = element_blank() ) + 
  ylab('Total reads assigned')
plot

ggave(plot, "comparison_plot.pdf", height = 5, width = 5)
