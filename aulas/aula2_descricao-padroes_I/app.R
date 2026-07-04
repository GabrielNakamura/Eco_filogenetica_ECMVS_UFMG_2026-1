library(shiny)
library(geiger)
library(ape)
library(phytools)

# UI
ui <- fluidPage(
  titlePanel("Evolução de atributos"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("n_lines", "Número de Espécies (Linhagens):", min = 5, max = 25, value = 10, step = 1),
      checkboxInput("show_mean", "Média do atributo no tempo presente", value = FALSE),
      hr(),
      sliderInput("rho", "Rho de Grafen (Distorção do Tempo):", min = 0.1, max = 3.0, value = 1.0, step = 0.1),
      p(tags$small("rho = 1: Padrão | rho < 1: Nós próximos à raiz | rho > 1: Nós próximos às pontas")),
      hr(),
      sliderInput("sig2_1", "Taxa de Evolução Painel 1 (sigma^2):", min = 0.01, max = 0.5, value = 0.02, step = 0.01),
      sliderInput("sig2_2", "Taxa de Evolução Painel 2 (sigma^2):", min = 0.01, max = 0.5, value = 0.10, step = 0.01),
      sliderInput("sig2_3", "Taxa de Evolução Painel 3 (sigma^2):", min = 0.01, max = 0.5, value = 0.40, step = 0.01),
      hr(),
      actionButton("go", "Simular")
    ),
    
    mainPanel(
      plotOutput("bmPlot", height = "1050px") 
    )
  )
)

# Server
server <- function(input, output, session) {
  
  sim_data <- eventReactive(c(input$go, input$n_lines, input$rho, input$sig2_1, input$sig2_2, input$sig2_3), {
    lines <- input$n_lines
    steps <- 100 
    
    # Gerar árvore inicial e aplicar distorção de Grafen
    raw_tree <- pbtree(n = lines, scale = 1, quiet = TRUE)
    raw_tree <- rotateNodes(raw_tree, "all") 
    
    distorted_tree <- compute.brlen(raw_tree, method = "Grafen", power = input$rho)
    distorted_tree$edge.length <- distorted_tree$edge.length / max(nodeHeights(distorted_tree))
    
    # Simulações Independentes Lineares (Linha 1)
    mat1 <- replicate(lines, c(0, cumsum(rnorm(steps, mean = 0, sd = sqrt(input$sig2_1)))))
    mat2 <- replicate(lines, c(0, cumsum(rnorm(steps, mean = 0, sd = sqrt(input$sig2_2)))))
    mat3 <- replicate(lines, c(0, cumsum(rnorm(steps, mean = 0, sd = sqrt(input$sig2_3)))))
    
    colnames(mat1) <- distorted_tree$tip.label
    colnames(mat2) <- distorted_tree$tip.label
    colnames(mat3) <- distorted_tree$tip.label
    
    # 3. Simulações dos atributos nas filogenias
    phy_sim1 <- fastBM(distorted_tree, sig2 = input$sig2_1, internal = TRUE)
    phy_sim2 <- fastBM(distorted_tree, sig2 = input$sig2_2, internal = TRUE)
    phy_sim3 <- fastBM(distorted_tree, sig2 = input$sig2_3, internal = TRUE)
    
    # Sincronização global para os gráficos lineares e mapas de árvores
    y_limits_global <- range(c(mat1, mat2, mat3, phy_sim1, phy_sim2, phy_sim3))
    x_limits_linear <- c(0, steps * 1.15)
    
    return(list(
      tree = distorted_tree, steps = steps, 
      x_limits_linear = x_limits_linear, y_limits_global = y_limits_global,
      mat1 = mat1, mat2 = mat2, mat3 = mat3,
      phy_sim1 = phy_sim1, phy_sim2 = phy_sim2, phy_sim3 = phy_sim3
    ))
  })
  
  output$bmPlot <- renderPlot({
    d <- sim_data()
    n_tips <- length(d$tree$tip.label)
    
    lineage_colors <- rainbow(n_tips, alpha = 0.8)
    names(lineage_colors) <- d$tree$tip.label
    
    # Grid de Layout: 3 Linhas x 3 Colunas 
    # Linha 1: 1, 2, 3 (Lineares)
    # Linha 2: 4, 5, 6 (Filogenias Clássicas)
    # Linha 3: 7, 8, 9 (Fenogramas Livres - Posicionados por último)
    layout(matrix(1:9, nrow = 3, byrow = TRUE))
    
    # --- LINHA 1: MODELO LINEAR INDEPENDENTE ---
    par(mar = c(4, 4.5, 3, 1.5))
    draw_linear_panel <- function(mat, sig2_val) {
      plot(0:d$steps, mat[, 1], type = "n", xlim = d$x_limits_linear, ylim = d$y_limits_global, 
           xlab = "Gerações", ylab = if(sig2_val == input$sig2_1) "Valor do Traço" else "", 
           main = paste("Linear (sig2 =", sig2_val, ")"), las = 1)
      grid(col = "lightgray", lty = "dotted")
      
      for(tip in d$tree$tip.label) {
        lines(0:d$steps, mat[, tip], col = lineage_colors[tip], lwd = 1.3)
        points(d$steps, mat[nrow(mat), tip], pch = 21, bg = lineage_colors[tip], col = "white", cex = 1.3)
        text(d$steps, mat[nrow(mat), tip], labels = tip, pos = 4, cex = 0.9, col = "black")
      }
      if (input$show_mean) abline(h = mean(mat[nrow(mat), ]), col = "black", lwd = 2.5, lty = "dashed")
    }
    draw_linear_panel(d$mat1, input$sig2_1)
    draw_linear_panel(d$mat2, input$sig2_2)
    draw_linear_panel(d$mat3, input$sig2_3)
    
    # --- LINHA 2: FILOGENIA COM distorção de GRAFEN ---
    par(mar = c(4, 2, 4, 3.5))
    draw_phylo_panel <- function(phy_sim, sig2_val) {
      tips_only <- phy_sim[d$tree$tip.label]
      cmap <- contMap(d$tree, tips_only, plot = FALSE)
      cmap <- setMap(cmap, c("darkgray", "lightgray", "darkgray")) 
      
      plot(cmap, type = "phylogram", direction = "rightwards", fsize = 1.0, lwd = 2, 
           main = paste("Árvore (sig2 =", sig2_val, ")"), legend = FALSE, mar = c(4, 2, 4, 3.5))
      
      tiplabels(pch = 21, bg = lineage_colors[d$tree$tip.label], col = "white", cex = 1.5)
      axisPhylo(side = 1)
      mtext("Tempo decorrido", side = 1, line = 2, cex = 0.8)
    }
    draw_phylo_panel(d$phy_sim1, input$sig2_1)
    draw_phylo_panel(d$phy_sim2, input$sig2_2)
    draw_phylo_panel(d$phy_sim3, input$sig2_3)
    
    # --- LINHA 3: FENOGRAMAS (REORGANIZADO PARA O FINAL COM EIXO Y INDEPENDENTE) ---
    par(mar = c(4.5, 4.5, 3, 1.5))
    draw_phenogram_panel <- function(phy_sim, sig2_val) {
      # Ajuste crítico: ylim agora captura o intervalo local de cada simulação específica,
      # aplicando um zoom automático independente para cada painel de fenograma.
      local_y_limits <- range(phy_sim)
      
      phenogram(d$tree, phy_sim, ylim = local_y_limits, 
                xlab = "Tempo decorrido", ylab = if(sig2_val == input$sig2_1) "Valor do atributo" else "", 
                main = paste("Fenograma Zoom (sig2 =", sig2_val, ")"), colors = "darkslateblue", las = 1, fsize = 0.9)
      
      # Plotagem das bolinhas sincronizadas de identificação no tempo final 1.0
      for(tip in d$tree$tip.label) {
        val_final <- phy_sim[tip]
        points(1.0, val_final, pch = 21, bg = lineage_colors[tip], col = "white", cex = 1.3)
      }
    }
    draw_phenogram_panel(d$phy_sim1, input$sig2_1)
    draw_phenogram_panel(d$phy_sim2, input$sig2_2)
    draw_phenogram_panel(d$phy_sim3, input$sig2_3)
  })
}

# Executa o app
shinyApp(ui = ui, server = server)
