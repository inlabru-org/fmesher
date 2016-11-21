## shinymesh.R
##
##   Copyright (C) 2016, Finn Lindgren
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

meshbuilder.app <- function() {

    meshbuilder.ui <- fluidPage(
        sidebarLayout(
            sidebarPanel(
                sliderInput(inputId="n", label="Seed points",
                            value=10, min=1, max=1000),
                actionButton(inputId="update.loc", label="Generate seed points"),
##                checkboxInput("use.boundary", label = "Boundary", value=FALSE),
                sliderInput(inputId="offset", label="offset",
                            value=c(0.1, 0.3), min=0.05, max=3),
                sliderInput(inputId="max.edge", label="max.edge",
                            value=c(0.05,0.5), min=0.01, max=2),
                sliderInput(inputId="min.angle", label="min.angle",
                            value=c(21,30), min=1, max=35),
                hr(),
                checkboxInput("assess", label = "Assessment", value=FALSE),
                radioButtons("assessment", label = "Type",
                             choices = list("Mesh SD" = "coarse",
                                            "SD ratio (mesh/fine)" = "rel",
                                            "Fine SD" = "fine"),
                             selected = "rel"),
                checkboxInput("sd.bound", label = "SD bound", value=FALSE),
                checkboxGroupInput("overlay", label="Overlay",
                                   choices = list("Mesh" = "mesh", "Fine" = "finemesh"),
                                   selected= "mesh"),
                sliderInput(inputId="corr.range", label="corr.range",
                            value=0.2, min=0.01, max=4),
                hr(),
                tableOutput("meta"),
                actionButton(inputId="plot.save", label="Save plot"),
                textOutput("plot.filename")
            ),
            mainPanel(
                plotOutput("meshplot", height="auto")
            )
        )
    )

    meshbuilder.server <- function(input, output, session) {
        plot.filename <- eventReactive(c(input$plot.save), {
            tempfile("meshbuilder_", getwd(), ".pdf")
        })
        output$plot.filename <- renderText({
            plot.filename()
        })

        observeEvent(c(input$plot.save), {
            if (input$plot.save > 0) { ## Don't save on app initialisation!
                pdf(plot.filename(), width=4, height=4)
                eval(plot.expr)
                dev.off()
                plot.filename()
            }
        })
        loc <- eventReactive(c(input$n, input$update.loc), {
            matrix(runif(input$n*2), input$n, 2)
        })
        boundary <- reactive({
            bnd <- list(
                INLA::inla.nonconvex.hull(loc(), input$offset[1]),
                INLA::inla.nonconvex.hull(loc(), input$offset[2]))
##            if (input$use.boundary) {
##                bnd[[2]] <- list(bnd[[2]],
##                                 inla.mesh.segment(matrix(c(0,1,1,0, 0,0,1,1), 4, 2),
##                                                   is.bnd=TRUE))
##            }
            bnd
        })
##        interior <- reactive({
##            int <- list(INLA::inla.mesh.segment(loc=matrix(c(0.2,0.8,0.5,0.5, 0.5,0.5,0.2,0.8), 4,2),
##                                                idx=matrix(c(1,3,2,4),2,2),
##                                                is.bnd=FALSE))
##            int
##        })
        mesh <- reactive({
            INLA::inla.mesh.2d(loc(),
                               boundary=boundary(),
##                               interior=interior(),
                               max.edge=input$max.edge,
                               min.angle=rev(input$min.angle),
                               cutoff=min(0.05, min(input$max.edge[1]/4, input$offset[1]/4)))
        })
        mesh.edgelengths <- reactive({
            m <- mesh()
            sqrt(c(Matrix::rowSums((m$loc[m$graph$tv[,2],] - m$loc[m$graph$tv[,1],])^2),
                   Matrix::rowSums((m$loc[m$graph$tv[,3],] - m$loc[m$graph$tv[,2],])^2),
                   Matrix::rowSums((m$loc[m$graph$tv[,1],] - m$loc[m$graph$tv[,3],])^2)))
        })
        finemesh.edgelengths <- reactive({
            m <- finemesh()
            if (is.null(m)) {
                c(NA)
            } else {
                sqrt(c(Matrix::rowSums((m$loc[m$graph$tv[,2],] - m$loc[m$graph$tv[,1],])^2),
                       Matrix::rowSums((m$loc[m$graph$tv[,3],] - m$loc[m$graph$tv[,2],])^2),
                       Matrix::rowSums((m$loc[m$graph$tv[,1],] - m$loc[m$graph$tv[,3],])^2)))
            }
        })
        finemesh <- reactive({
            if (!input$assess) {
                NULL
            } else {
                INLA::inla.mesh.2d(loc(),
                                   boundary=list(INLA::inla.mesh.boundary(mesh())),
##                                   interior=interior(),
                                   max.edge=min(mesh.edgelengths())/4,
                                   min.angle=21,
                                   max.n=mesh()$n*16,
                                   cutoff=0)
            }
        })
        mesh.proj <- reactive({
            INLA::inla.mesh.projector(mesh(), dims=c(500, 500))
        })
        finemesh.proj <- reactive({
            if (!input$assess) {
                NULL
            } else {
                INLA::inla.mesh.projector(finemesh(), lattice=mesh.proj()$lattice)
            }
        })
        mesh.S <- reactive({
            if (!input$assess) {
                NULL
            } else {
                spde <- INLA::inla.spde2.pcmatern(mesh(),
                                                  prior.range=c(input$corr.range, 0.5),
                                                  prior.sigma=c(1, 0.5))
                Q <- INLA::inla.spde.precision(spde, theta=log(c(input$corr.range, 1)))
                INLA::inla.qinv(Q, reordering=INLA::inla.reorderings())
            }
        })
        finemesh.S <- reactive({
            if (!input$assess) {
                NULL
            } else {
                spde <- INLA::inla.spde2.pcmatern(finemesh(),
                                                  prior.range=c(input$corr.range, 0.5),
                                                  prior.sigma=c(1, 0.5))
                Q <- INLA::inla.spde.precision(spde, theta=log(c(input$corr.range, 1)))
                INLA::inla.qinv(Q, reordering=INLA::inla.reorderings())
            }
        })
        mesh.sd <- reactive({
            if (!input$assess) {
                NULL
            } else {
                proj <- mesh.proj()
                if (input$sd.bound) {
                    INLA::inla.mesh.project(proj, field=diag(mesh.S())^0.5)
                } else {
                    v <- Matrix::rowSums(proj$proj$A * (proj$proj$A %*% mesh.S()))
                    v[!proj$proj$ok] <- NA
                    matrix(v^0.5, length(proj$x), length(proj$y))
                }
            }
        })
        finemesh.sd <- reactive({
            if (!input$assess) {
                NULL
            } else {
                proj <- finemesh.proj()
                if (input$sd.bound) {
                    INLA::inla.mesh.project(proj, field=diag(finemesh.S())^0.5)
                } else {
                    v <- Matrix::rowSums(proj$proj$A * (proj$proj$A %*% finemesh.S()))
                    v[!proj$proj$ok] <- NA
                    matrix(v^0.5, length(proj$x), length(proj$y))
                }
            }
        })
        
        meta <- reactive({
            m <- rbind(c("Mesh", mesh()$n, range(mesh.edgelengths())))
            m2 <- finemesh()
            if (is.null(m2)) {
                m <- rbind(m, c("Fine", NA, NA, NA))
            } else {
                m <- rbind(m, c("Fine", finemesh()$n, range(finemesh.edgelengths())))
            }
            colnames(m) <- c("Resolution", "nV", "minE", "maxE")
            as.data.frame(m)
        })
        
        plot.expr <- quote({
            if (!input$assess) {
                plot(mesh(), asp=1, main="")
            } else {
                col <- colorRampPalette(c("blue", "white", "red"))
                n.col <- 1+16
                zlim.abs <- range(c(range(mesh.sd(), na.rm=TRUE),
                                    range(finemesh.sd(), na.rm=TRUE)))
                zlim.abs <- exp(max(abs(log(zlim.abs))))
                zlim.abs <- 1 + c(-1,1)*(zlim.abs-1)
                zlim.rel <- range(mesh.sd()/finemesh.sd(), na.rm=TRUE)
                zlim.rel <- exp(max(abs(log(zlim.rel))))
                zlim.rel <- 1 + c(-1,1)*max(0.1, zlim.rel-1)
                if (input$assessment == "coarse") {
                    fields::image.plot(mesh.proj()$x, mesh.proj()$y, (mesh.sd()),
                                       zlim=zlim.abs,
                                       col=col(n.col), asp=1, main="")
                } else if (input$assessment == "rel") {
                    fields::image.plot(mesh.proj()$x, mesh.proj()$y, (mesh.sd()/finemesh.sd()),
                                       zlim=zlim.rel,
                                       col=col(n.col), asp=1, main="")
                } else if (input$assessment == "fine") {
                    fields::image.plot(mesh.proj()$x, mesh.proj()$y, (finemesh.sd()),
                                       zlim=zlim.abs,
                                       col=col(n.col), asp=1, main="")
                }
                if ("mesh" %in% input$overlay) {
                    plot(mesh(), add=TRUE)
                }
                if ("finemesh" %in% input$overlay) {
                    plot(finemesh(), add=TRUE)
                }
            }
        })

        output$meshplot <-
            renderPlot(plot.expr, quoted=TRUE, height=1000, units="px")

        output$meta <- renderTable({meta()})
    }

    shinyApp(ui=meshbuilder.ui, server=meshbuilder.server)
}


##' @import shiny
##' @import Matrix
##' @export
meshbuilder <- function() {
    requireNamespace("INLA")
    requireNamespace("fields")
    library("INLA") ## Workaround for INLA bug prior to 2016-11-15

    runApp(meshbuilder.app())
}
