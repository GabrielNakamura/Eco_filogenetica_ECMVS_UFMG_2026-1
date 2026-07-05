library(shiny)
library(phytools)
library(nlme) # Necessário para rodar a função gls() e corBrownian()

# Define UI
ui <- fluidPage(
  titlePanel("Evolução de atributos, Movimento Browniano, PIC e PGLS"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("n_lines", "Número de Espécies (Linhagens):", min = 6, max = 25, value = 12, step = 1),
      checkboxInput("show_mean", "Mostrar Média Geral", value = FALSE),
      hr(),
      sliderInput("rho", "Distorção do Tempo:", min = 0.1, max = 3.0, value = 1.0, step = 0.1),
      p(tags$small("rho = 1: Padrão | rho < 1: Nós próximos à raiz | rho > 1: Nós próximos às pontas")),
      hr(),
      sliderInput("sig2_1", "Taxa de Evolução Painel 1 (sigma^2):", min = 0.01, max = 0.5, value = 0.02, step = 0.01),
      sliderInput("sig2_2", "Taxa de Evolução Painel 2 (sigma^2):", min = 0.01, max = 0.5, value = 0.10, step = 0.01),
      sliderInput("sig2_3", "Taxa de Evolução Painel 3 (sigma^2):", min = 0.01, max = 0.5, value = 0.40, step = 0.01),
      hr(),
      actionButton("go", "Simular")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Modelos Evolutivos", 
                 plotOutput("bmPlot", height = "1050px")
        ),
        tabPanel("Relação de Atributos & PIC", 
                 plotOutput("picPlot", height = "500px")
        ),
        tabPanel("Relação de Atributos & PGLS", 
                 plotOutput("traitPlot", height = "950px")
        )
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Bloco Reativo Principal de Simulação
  sim_data <- eventReactive(c(input$go, input$n_lines, input$rho, input$sig2_1, input$sig2_2, input$sig2_3), {
    lines <- input$n_lines
    steps <- 100 
    
    # 1. Gerar árvore e aplicar distorção de Grafen via compute.brlen (pacote ape)
    raw_tree <- pbtree(n = lines, scale = 1, quiet = TRUE)
    raw_tree <- rotateNodes(raw_tree, "all") 
    
    distorted_tree <- compute.brlen(raw_tree, method = "Grafen", power = input$rho)
    distorted_tree$edge.length <- distorted_tree$edge.length / max(nodeHeights(distorted_tree))
    
    # 2. Simulações Independentes Lineares (Aba 1)
    mat1 <- replicate(lines, c(0, cumsum(rnorm(steps, mean = 0, sd = sqrt(input$sig2_1)))))
    mat2 <- replicate(lines, c(0, cumsum(rnorm(steps, mean = 0, sd = sqrt(input$sig2_2)))))
    mat3 <- replicate(lines, c(0, cumsum(rnorm(steps, mean = 0, sd = sqrt(input$sig2_3)))))
    
    colnames(mat1) <- distorted_tree$tip.label
    colnames(mat2) <- distorted_tree$tip.label
    colnames(mat3) <- distorted_tree$tip.label
    
    # 3. Simulações Filogenéticas (Atributo A)
    phy_sim1_A <- fastBM(distorted_tree, sig2 = input$sig2_1, internal = TRUE)
    phy_sim2_A <- fastBM(distorted_tree, sig2 = input$sig2_2, internal = TRUE)
    phy_sim3_A <- fastBM(distorted_tree, sig2 = input$sig2_3, internal = TRUE)
    
    # 4. Simulação de um SEGUNDO Atributo Independente (Atributo B)
    phy_sim1_B <- fastBM(distorted_tree, sig2 = input$sig2_1, internal = TRUE)
    phy_sim2_B <- fastBM(distorted_tree, sig2 = input$sig2_2, internal = TRUE)
    phy_sim3_B <- fastBM(distorted_tree, sig2 = input$sig2_3, internal = TRUE)
    
    y_limits_global <- range(c(mat1, mat2, max(mat3), phy_sim1_A, phy_sim2_A, phy_sim3_A))
    x_limits_linear <- c(0, steps * 1.15)
    
    return(list(
      tree = distorted_tree, steps = steps, 
      x_limits_linear = x_limits_linear, y_limits_global = y_limits_global,
      mat1 = mat1, mat2 = mat2, mat3 = mat3,
      phy_sim1_A = phy_sim1_A, phy_sim2_A = phy_sim2_A, phy_sim3_A = phy_sim3_A,
      phy_sim1_B = phy_sim1_B, phy_sim2_B = phy_sim2_B, phy_sim3_B = phy_sim3_B
    ))
  })
  
  # --- RENDERIZAÇÃO ABA 1: PAINÉIS EVOLUTIVOS  ---
  output$bmPlot <- renderPlot({
    d <- sim_data()
    n_tips <- length(d$tree$tip.label)
    lineage_colors <- rainbow(n_tips, alpha = 0.8)
    names(lineage_colors) <- d$tree$tip.label
    
    layout(matrix(1:9, nrow = 3, byrow = TRUE))
    
    # Linha 1: Modelos Lineares
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
    
    # Linha 2: Árvores Filogenéticas Clássicas
    par(mar = c(4, 2, 4, 3.5))
    draw_phylo_panel <- function(phy_sim, sig2_val) {
      tips_only <- phy_sim[d$tree$tip.label]
      cmap <- contMap(d$tree, tips_only, plot = FALSE)
      cmap <- setMap(cmap, c("darkgray", "lightgray", "darkgray")) 
      plot(cmap, type = "phylogram", direction = "rightwards", fsize = 1.0, lwd = 2, 
           main = paste("Árvore (sig2 =", sig2_val, ")"), legend = FALSE, mar = c(4, 2, 4, 3.5))
      tiplabels(pch = 21, bg = lineage_colors[d$tree$tip.label], col = "white", cex = 1.5)
      axisPhylo(side = 1)
    }
    draw_phylo_panel(d$phy_sim1_A, input$sig2_1)
    draw_phylo_panel(d$phy_sim2_A, input$sig2_2)
    draw_phylo_panel(d$phy_sim3_A, input$sig2_3)
    
    # Linha 3: Fenogramas
    par(mar = c(4.5, 4.5, 3, 1.5))
    draw_phenogram_panel <- function(phy_sim, sig2_val) {
      local_y_limits <- range(phy_sim)
      phenogram(d$tree, phy_sim, ylim = local_y_limits, 
                xlab = "Tempo Histórico", ylab = if(sig2_val == input$sig2_1) "Valor do Traço" else "", 
                main = paste("Fenograma Zoom (sig2 =", sig2_val, ")"), colors = "darkslateblue", las = 1, fsize = 0.9)
      for(tip in d$tree$tip.label) {
        points(1.0, phy_sim[tip], pch = 21, bg = lineage_colors[tip], col = "white", cex = 1.3)
      }
    }
    draw_phenogram_panel(d$phy_sim1_A, input$sig2_1)
    draw_phenogram_panel(d$phy_sim2_A, input$sig2_2)
    draw_phenogram_panel(d$phy_sim3_A, input$sig2_3)
  })
  
  # --- RENDERIZAÇÃO ABA 2: CONTRASTES FILOGENÉTICOS INDEPENDENTES (PIC) ---
  output$picPlot <- renderPlot({
    d <- sim_data()
    
    # Configura o layout: 1 linha por 3 colunas para os três sigmas
    layout(matrix(1:3, nrow = 1))
    par(mar = c(5, 5, 4, 1.5))
    
    draw_pic_panel <- function(phy_A, phy_B, sig2_val) {
      x_vals <- phy_A[d$tree$tip.label]
      y_vals <- phy_B[d$tree$tip.label]
      
      # Calcular os contrastes independentes (PIC) utilizando a função nativa pic() do pacote ape
      pic_X <- pic(x_vals, d$tree)
      pic_Y <- pic(y_vals, d$tree)
      
      # Ajustar modelo linear passando OBRIGATORIAMENTE pela origem (Y ~ X - 1 ou Y ~ X + 0)
      pic_fit <- lm(pic_Y ~ pic_X + 0)
      pic_summary <- summary(pic_fit)
      
      beta_pic <- round(pic_summary$coefficients[1, 1], 4)
      p_pic <- round(pic_summary$coefficients[1, 4], 4)
      
      # Gerar gráfico de dispersão dos contrastes
      plot(pic_X, pic_Y, xlab = "Contrastes do Atributo X", ylab = "Contrastes do Atributo Y",
           main = paste("Regressão via PIC\nsigma^2 =", sig2_val), las = 1,
           pch = 21, bg = "forestgreen", col = "darkgreen", cex = 1.5)
      grid(col = "lightgray", lty = "dotted")
      
      # Força a linha de regressão a cruzar o ponto central (0,0)
      abline(pic_fit, col = "darkgreen", lwd = 2.5)
      abline(h = 0, v = 0, col = "gray", lty = "dashed", lwd = 0.8) # Linhas de referência da origem
      
      # Anotação com os parâmetros estatísticos do PIC
      mtext(paste0("PIC (Regressão na Origem):\nBeta = ", beta_pic, "\nP-valor = ", p_pic), 
            side = 3, line = -3, cex = 0.9, col = "darkgreen", font = 2, adj = 0.05)
    }
    
    draw_pic_panel(d$phy_sim1_A, d$phy_sim1_B, input$sig2_1)
    draw_pic_panel(d$phy_sim2_A, d$phy_sim2_B, input$sig2_2)
    draw_pic_panel(d$phy_sim3_A, d$phy_sim3_B, input$sig2_3)
  })
  
  # --- RENDERIZAÇÃO ABA 3: GLM CLÁSSICO vs PGLS (FILOGENÉTICO) ---
  output$traitPlot <- renderPlot({
    d <- sim_data()
    n_tips <- length(d$tree$tip.label)
    lineage_colors <- rainbow(n_tips, alpha = 0.8)
    names(lineage_colors) <- d$tree$tip.label
    
    layout(matrix(1:6, nrow = 2, byrow = TRUE))
    
    # Sub-rotina para o GLM clássico
    draw_glm_panel <- function(phy_A, phy_B, sig2_val) {
      par(mar = c(4.5, 5, 4, 1.5))
      x_vals <- phy_A[d$tree$tip.label]
      y_vals <- phy_B[d$tree$tip.label]
      
      df_data <- data.frame(X = x_vals, Y = y_vals)
      glm_fit <- glm(Y ~ X, data = df_data, family = gaussian)
      glm_summary <- summary(glm_fit)
      
      beta_val <- round(glm_summary$coefficients[2, 1], 4) 
      p_val <- round(glm_summary$coefficients[2, 4], 4)    
      
      lambda_X <- round(phylosig(d$tree, x_vals, method = "lambda")$lambda, 3)
      lambda_Y <- round(phylosig(d$tree, y_vals, method = "lambda")$lambda, 3)
      
      plot(x_vals, y_vals, type = "n", xlab = "Atributo X", ylab = "Atributo Y",
           main = paste("GLM Tradicional (Normal)\nsigma^2 =", sig2_val), las = 1)
      grid(col = "lightgray", lty = "dotted")
      abline(glm_fit, col = "darkred", lwd = 2.5, lty = "dashed")
      
      for(tip in d$tree$tip.label) {
        points(x_vals[tip], y_vals[tip], pch = 21, bg = lineage_colors[tip], col = "white", cex = 1.6)
        text(x_vals[tip], y_vals[tip], labels = tip, pos = 4, cex = 0.85, col = "darkgray")
      }
      
      mtext(paste0("GLM: Beta = ", 
                   beta_val, 
                   " | P-valor = ", 
                   p_val, 
                   "\nSinal Filo: Lambda X = ", 
                   lambda_X, 
                   " | Lambda Y = ",
                   lambda_Y), 
            side = 3, 
            line = -3,
            cex = 0.8, 
            col = "darkred",
            font = 2, 
            adj = 0.02)
    }
    
    # Sub-rotina para o PGLS filogenético
    draw_pgls_panel <- 
      function(phy_A, phy_B, sig2_val) {
      par(mar = c(4.5, 5, 4, 1.5))
        x_vals <- phy_A[d$tree$tip.label]
        y_vals <- phy_B[d$tree$tip.label]
        df_data <- data.frame(X = x_vals, Y = y_vals)
        rownames(df_data) <- d$tree$tip.label
        bm_cor <- corBrownian(phy = d$tree)
        pgls_fit <- gls(Y ~ X, data = df_data, correlation = bm_cor, method = "ML")
        pgls_summary <- summary(pgls_fit)
        beta_pgls <- round(pgls_summary$coefficients[2], 4)
        p_pgls <- round(pgls_summary$tTable[2, 4], 4)
        plot(x_vals, y_vals, 
             type = "n", xlab = "Atributo X",
             ylab = "Atributo Y", 
             main = paste("PGLS (corBrownian)\nsigma^2 =", sig2_val),
             las = 1)
        grid(col = "lightgray", lty = "dotted")
        abline(a = pgls_summary$coefficients[1], b = pgls_summary$coefficients[2], col = "darkblue", lwd = 2.5)
        for(tip in d$tree$tip.label) {
          points(x_vals[tip], y_vals[tip], pch = 21, bg = lineage_colors[tip], col = "white", cex = 1.6)
          text(x_vals[tip], y_vals[tip], labels = tip, pos = 4, cex = 0.85, col = "darkgray")
        }
        mtext(paste0("PGLS (BM): Beta = ", beta_pgls, " | P-valor = ", p_pgls), 
              side = 3, line = -2, cex = 0.85, col = "darkblue", font = 2, adj = 0.02)
      }
    # Executa as sublinhas do PGLS
    draw_glm_panel(d$phy_sim1_A, d$phy_sim1_B, input$sig2_1)
    draw_glm_panel(d$phy_sim2_A, d$phy_sim2_B, input$sig2_2)
    draw_glm_panel(d$phy_sim3_A, d$phy_sim3_B, input$sig2_3)
    draw_pgls_panel(d$phy_sim1_A, d$phy_sim1_B, input$sig2_1)
    draw_pgls_panel(d$phy_sim2_A, d$phy_sim2_B, input$sig2_2)
    draw_pgls_panel(d$phy_sim3_A, d$phy_sim3_B, input$sig2_3)
  }
  )
  }

shinyApp(ui = ui, server = server)      