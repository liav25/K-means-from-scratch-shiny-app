#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(dplyr)
library(plotly)
library(shiny)
library(shinydashboard)

med_dat = read.delim ("gtex.gct", skip = 2, row.names=c(1), header = TRUE)
gen_names = med_dat[,1]
med_dat = med_dat[,-c(1)]

tissues = colnames(med_dat)

med_dat_t = t(med_dat)

med_dat_t_clean = as.data.frame(med_dat_t[,which(apply(med_dat_t, 2, sd)!=0)])

#create pc analysis
pca = prcomp(med_dat_t_clean, center = TRUE ,scale = TRUE)

#check how much columns we need:

#In this part we checked how many dimension should we use for pca.
#We can do so by checking the sum of explained variance by the Principle Components. 

pca.var = pca$sdev^2
pca.var.per = round(pca.var/sum(pca.var)*100, 3)
pca.cumsum = cumsum(pca.var.per)
df = pca.var.per %>% cbind(cumsum=pca.cumsum) %>%  `colnames<-`(c("PER", "CUMU")) %>% as.data.frame()
ggplot(df) + geom_col(aes(x=1:53, y=PER, fill="Explained Variance")) + geom_point(aes(x=1:53, y=CUMU,color='Cumulative Sum of Precentage') ) +
  theme_light() + labs(title="Cumulative Scree Plot", x="", y="Precentage") + 
  scale_fill_manual(values=c("#00AFBB"), name="") + scale_color_manual(name="", values=c("orangered"))

# We found that for 40 Principle component we get 98.6% of the variance

med_dat_40d = as.data.frame(pca$x[,1:40])

#create 3-dimensional dataframe
med_dat_3d = as.data.frame(pca$x[,1:3])
#create 2-dimensional dataframe
med_dat_2d = as.data.frame(pca$x[,1:2])


# this function initials the first centroids randomally from our data
# param k: number of centroids
# param data: our data
# return: a vector of centroids
init_centroids = function(k, data) {
  cent = data[sample(1:nrow(data), k), ]
  return(cent)
}

# this function assign split data to clusters by cluster centroids list
# param centroids: a list of centroids
# param data: dataframe
# return: a list of clusters (by index)
assign_to_clusters = function(centroids, data) {
  colnames(centroids) = colnames(data)
  distance = dist(rbind(data,centroids))
  n_dat = nrow(data)
  n_cen = nrow(centroids)
  d = as.matrix(distance)[(n_dat + 1):(n_dat+n_cen),1:n_dat]
  results = apply(d,2,function(x) which(x == min(x)))
  clusters = as.numeric(unlist(results))  
  return(clusters)
}


# this function recalculate centroids of clusters by calculating the new mean of cluster
# param data: the data
# param clus_lst: the clusters list (indices of clusters)
# param k: number of clusters
# return: new centroids of current clusters
recalc_centroids = function(data, clus_lst, k) {
  d = data.frame(matrix(0, ncol = dim(data)[2],nrow=k)) 
  for (i in 1:k) {
    cent = colMeans(data[clus_lst==i,])
    d[i,] = cent
  }
  
  return(d)
}

check_diff = function(cent_after, cent_before, eps){
  dif = as.matrix(cent_after - cent_before)
  score = apply(dif, 1, norm, type="2")
  if (all(score<=eps)){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

# this is a k means algorithm from scratch
k_means = function(k, data, eps=1e-3, max_iter=10) {
  #1
  #2+3
  cent0 = init_centroids(k, data)
  
  
  clusters0 = assign_to_clusters(centroids = cent0, data = data)
  #4
  cent = recalc_centroids(data=data,clus_lst = clusters0 , k=k)
  #5
  is_diff = check_diff(cent_before = cent0, cent_after = cent, eps=eps)
  iter = 0
  while( (is_diff) && (iter <= max_iter) ){
    iter = iter + 1
    clusters = assign_to_clusters(centroids = cent, data=data)
    cent0 = cent
    cent = recalc_centroids(data=data,clus_lst = clusters , k=k)
    clusters = assign_to_clusters(centroids = cent, data=data)
    is_diff = check_diff(cent_before = cent0, cent_after = cent, eps=eps)
  }
  return(clusters)
}


############
############
############
############
############
############


# Define UI for random distribution app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("KMeans Algorithm Comparison - 2D and 3D Genes Data"),
  
  h4("asdadsad<br>asdasdasdasd"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # br() element to introduce extra vertical spacing ----
      br(),
      
      # Input: Slider for the number of observations to generate ----
      sliderInput("slider",
                  "Number of clusters:",
                  value = 4,
                  min = 2,
                  max = 10)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("2D Plot", plotlyOutput("plot")),
                  tabPanel("3D Plot", plotlyOutput("summary"))
                  
      )
      
    )
  )
)
# Define server logic for random distribution app ----
server <- function(input, output) {
  
  output$plot <- renderPlotly({
    clusters = k_means(k=input$slider, data=med_dat_40d)
    results_df = cbind(med_dat_2d,as.factor(clusters))
    colnames(results_df) = c("PC1","PC2","color")
    
    plot_ly(data=results_df, x=~PC1, y=~PC2,  
            type="scatter", mode="markers", color=~color, text=tissues)     
  })
  
  output$summary <- renderPlotly({
    clusters = k_means(k=input$slider, data=med_dat_40d)
    
    results_df = cbind(med_dat_3d,as.factor(clusters))
    colnames(results_df) = c("PC1","PC2","PC3","color")
    
    plot_ly(data=results_df, x=~PC1, y=~PC2, z=~PC3, 
               type="scatter3d", mode="markers", color=~color, text=tissues) 
  })
  
}



shinyApp(ui = ui, server = server)

