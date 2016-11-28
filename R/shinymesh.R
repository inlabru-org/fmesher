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
                h4("Settings for inla.nonconvex.hull()"),
                sliderInput(inputId="offset", label="offsets for automatic boundaries",
                            value=c(0.1, 0.3), min=0.05, max=2),
                h4("Settings for inla.mesh.2d()"),
                sliderInput(inputId="max.edge", label="max.edge",
                            value=c(0.05,0.5), min=0.01, max=2),
                sliderInput(inputId="cutoff", label="cutoff",
                            value=0.01, min=0, max=0.2),
                sliderInput(inputId="min.angle", label="min.angle",
                            value=c(21,30), min=1, max=35),
                hr(),
                h4("Mesh resolution assessment"),
                fluidRow(
                    column(6,
                           checkboxInput("assess", label = "Active", value=FALSE),
                           radioButtons("assessment", label = "Type",
                                        choices = list("Mesh SD" = "coarse",
                                                       "SD ratio (mesh/fine)" = "rel",
                                                       "Fine SD" = "fine"),
                                        selected = "rel")
                           ),
                    column(6,
                           checkboxInput("sd.bound", label = "SD bound", value=FALSE),
                           checkboxGroupInput("overlay", label="Overlay",
                                              choices = list("Mesh" = "mesh",
                                                             "Fine" = "finemesh",
                                                             "Correlations" ="corr"),
                                              selected= "mesh")
                           )
                ),
                sliderInput(inputId="corr.range", label="Spatial correlation range",
                            value=0.2, min=0.01, max=4),
                sliderInput(inputId="corr.nu", label="Matern smoothness (nu)",
                            value=1, min=0.01, max=1),
                hr(),
                actionButton(inputId="plot.save", label="Save plot"),
                verbatimTextOutput("plot.filename"),
                width=3
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Input",
                             fluidRow(
                                 column(4,
                                        h4("Random seed points"),
                                        checkboxGroupInput(inputId="loc.usage",
                                                     label="Use random seed points for",
                                                     choices=list("Domain boundary definitions"="boundary",
                                                                  "Mesh seed vertices"="mesh"),
                                                     selected="boundary"),
                                        sliderInput(inputId="input.loc.n", label="Seed points",
                                                    value=10, min=1, max=1000),
                                        actionButton(inputId="update.loc",
                                                     label="Generate seed points"),
                                        h4("User defined points"),
                                        selectInput("boundary.loc.name",
                                                    label="Boundary domain points variable name(s)",
                                                    multiple=TRUE,
                                                    choices=with(globalenv(), ls())),
                                        selectInput("mesh.loc.name",
                                                    label="Mesh vertex seed points variable name(s)",
                                                    multiple=TRUE,
                                                    choices=with(globalenv(), ls())),
                                        h4("User defined boundaries"),
                                        selectInput("boundary1.name",
                                                    label="Inner boundary variable name(s)",
                                                    multiple=TRUE,
                                                    choices=with(globalenv(), ls())),
                                        selectInput("boundary2.name",
                                                    label="Outer boundary variable name(s)",
                                                    multiple=TRUE,
                                                    choices=with(globalenv(), ls())),
                                        verbatimTextOutput("crs.strings"),
                                        verbatimTextOutput("input.coordinates")
                                        ),
                                 column(8,
                                        plotOutput("inputplot", height="auto",
                                                   hover="inputplot.hover")
                                        )
                                 )
                             ),
                    tabPanel("Display",
                             plotOutput("meshplot", height="auto",
                                        click="meshplot.click"
                                        ),
                             fluidRow(
                                 column(4,
                                        tableOutput("meta")
                                        ),
                                 column(4,
                                        tableOutput("SD.information")
                                        ),
                                 column(4,
                                        tableOutput("corr.information")
                                        )
                             ),
                             verbatimTextOutput("loc.click.coordinates")
                             ),
                    tabPanel("Code",
                             verbatimTextOutput("code")
                             ),
                    tabPanel("Debug",
                             verbatimTextOutput("debug")
                             ),
                    id="panel",
                    selected="Input"
                )
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
                eval(meshplot.expr)
                dev.off()
                plot.filename()
            }
        })
        loc.seed <- eventReactive(c(input$update.loc), {
            as.integer(runif(1, min=0, max=.Machine$integer.max))
        })
        input.loc.crs <- reactive({
            sp <- c(boundary.loc.input(), mesh.loc.input())
            if (!is.null(sp) && (length(sp) > 0)) {
                stopifnot(identicalCRS(sp))
                CRS(proj4string(sp[[1]]))
            } else {
                CRS()
            }
        })
        loc <- reactive({
            set.seed(loc.seed())
            out <- SpatialPoints(matrix(runif(2*input$input.loc.n) *
                                        c(diff(userinput.xlim()), diff(userinput.ylim())) +
                                        c(userinput.xlim()[1], userinput.ylim()[1]),
                                        input$input.loc.n, 2, byrow=TRUE),
                                 input.loc.crs())
            attr(out, "code") <- paste0(
                "set.seed(", loc.seed(), ")\n",
                "loc <- SpatialPoints(matrix(runif(2*", input$input.loc.n, "), ",
                input$input.loc.n, ", 2, byrow=TRUE))\n"
            )
            out
        })
        boundary.loc.input <- reactive({
            out <- lapply(input$boundary.loc.name,
                          function(name) {
                              if (length(intersect(class(globalenv()[[name]]),
                                                   c("SpatialPoints", "SpatialPointsDataFrame",
                                                     "SpatialLines", "SpatialLinesDataFrame",
                                                     "SpatialPolygons", "SpatialPolygonsDataFrame"))) > 0) {
                                  coord <- coordinates(globalenv()[[name]])
                                  crs <- CRS(proj4string(globalenv()[[name]]))
                              } else {
                                  coord <- as.matrix(globalenv()[[name]])
                                  crs <- CRS()
                              }
                              SpatialPoints(coord, crs)
                          })
            names(out) <- input$boundary.loc.name
            if (length(out) > 0) {
                out
            } else {
                NULL
            }
        })
        mesh.loc.input <- reactive({
            out <- lapply(input$mesh.loc.name,
                          function(name) {
                              if (length(intersect(class(globalenv()[[name]]),
                                                   c("SpatialPoints", "SpatialPointsDataFrame",
                                                     "SpatialLines", "SpatialLinesDataFrame",
                                                     "SpatialPolygons", "SpatialPolygonsDataFrame"))) > 0) {
                                  coord <- coordinates(globalenv()[[name]])
                                  crs <- CRS(proj4string(globalenv()[[name]]))
                              } else {
                                  coord <- as.matrix(globalenv()[[name]])
                                  crs <- CRS()
                              }
                              SpatialPoints(coord, crs)
                          })
            names(out) <- input$mesh.loc.name
            if (length(out) > 0) {
                out
            } else {
                NULL
            }
        })

        boundary.loc <- reactive({
            out <- NULL
            code <- ""
            sp <- boundary.loc.input()
            num <- length(sp) + ("boundary" %in% input$loc.usage)
            if (num > 0) {
                code <- "boundary.loc <- "
                if (num == 1) {
                    glue <- ""
                } else {
                    glue <- "rbind("
                }
                if (length(sp) == 1) {
                    name <- input$boundary.loc.name[1]
                    out <- sp[[name]]
                    code <- paste0(code, glue, name)
                    glue <- ", "
                } else  if (length(sp) > 1) {
                    out <- do.call(rbind, lapply(input$boundary.loc.name,
                                                 function(name) sp[[name]]))
                    code <- paste0(code, glue, paste0(input$boundary.loc.name, collapse=", "))
                    glue <- ", "
                }
                if ("boundary" %in% input$loc.usage) {
                    if (is.null(out)) {
                        out <- loc()
                    } else {
                        out <- rbind(out, loc())
                    }
                    code <- paste0(code, glue, "loc")
                }
                if (num > 1) {
                    code <- paste0(code, ")")
                }
                attr(out, "code") <- paste0(code, "\n")
            }
            out
        })
        mesh.loc <- reactive({
            out <- NULL
            code <- ""
            sp <- mesh.loc.input()
            num <- length(sp) + ("mesh" %in% input$loc.usage)
            if (num > 0) {
                code <- "mesh.loc <- "
                if (num == 1) {
                    glue <- ""
                } else {
                    glue <- "rbind("
                }
                if (length(sp) == 1) {
                    name <- input$mesh.loc.name[1]
                    out <- sp[[name]]
                    code <- paste0(code, glue, name)
                    glue <- ", "
                } else  if (length(sp) > 1) {
                    out <- do.call(rbind, lapply(input$mesh.loc.name,
                                                 function(name) sp[[name]]))
                    code <- paste0(code, glue, paste0(input$mesh.loc.name, collapse=", "))
                    glue <- ", "
                }
                if ("mesh" %in% input$loc.usage) {
                    if (is.null(out)) {
                        out <- loc()
                    } else {
                        out <- rbind(out, loc())
                    }
                    code <- paste0(code, glue, "loc")
                }
                if (num > 1) {
                    code <- paste0(code, ")")
                }
                attr(out, "code") <- paste0(code, "\n")
            }
            out
        })
        boundary <- reactive({
            loc1 <- boundary.loc()
            if (is.null(loc1)) {
                bnd <- NULL
            } else {
                bnd <- list(
                    INLA::inla.nonconvex.hull(coordinates(loc1), input$offset[1]),
                    INLA::inla.nonconvex.hull(coordinates(loc1), input$offset[2]))
                attr(bnd, "code") <- paste0(
                    "boundary <- list(inla.nonconvex.hull(coordinates(boundary.loc), ", input$offset[1], "),\n",
                    "                 inla.nonconvex.hull(coordinates(boundary.loc), ", input$offset[2], "))\n"
                )
            }
            bnd
        })
        mesh <- reactive({
            if (input$panel == "Input") {
                return(NULL)
            }
            message("mesh")
            loc1 <- mesh.loc()
            bnd1 <- boundary()
            if (is.null(loc1)) {
                loc.code <- ""
            } else {
                loc.code <- "loc=mesh.loc,
                     "
            }
            if (is.null(bnd1)) {
                bnd.code <- ""
            } else {
                bnd.code <- "boundary=boundary,
                     "
            }
            if (is.null(loc1) && is.null(bnd1)) {
                out <- NULL
            } else {
                out <- INLA::inla.mesh.2d(loc=loc1,
                                          boundary=bnd1,
                                          max.edge=input$max.edge,
                                          min.angle=rev(input$min.angle),
                                          max.n=c(48000, 16000),
                                          cutoff=input$cutoff)
                attr(out, "code") <- paste0(
"mesh <- inla.mesh.2d(", loc.code, bnd.code,
                    "max.edge=c(", input$max.edge[1], ", ", input$max.edge[2], "),
                     min.angle=c(", input$min.angle[2], ", ", input$min.angle[1], "),
                     max.n=c(48000, 16000)
                     cutoff=", input$cutoff, ")\n")
            }
            out
        })
        mesh.edgelengths <- reactive({
            m <- mesh()
            if (is.null(m)) {
                c(NA)
            } else {
                sqrt(c(Matrix::rowSums((m$loc[m$graph$tv[,2],] - m$loc[m$graph$tv[,1],])^2),
                       Matrix::rowSums((m$loc[m$graph$tv[,3],] - m$loc[m$graph$tv[,2],])^2),
                       Matrix::rowSums((m$loc[m$graph$tv[,1],] - m$loc[m$graph$tv[,3],])^2)))
            }
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
            if (!input$assess || is.null(mesh())) {
                out <- NULL
            } else {
                message("finemesh")
                out <- INLA::inla.mesh.2d(loc=mesh()$loc,
                                          boundary=list(INLA::inla.mesh.boundary(mesh())),
                                          max.edge=min(mesh.edgelengths())/2,
                                          min.angle=21,
                                          max.n=64000,
                                          cutoff=0)
                attr(out, "code") <- paste0(
"finemesh <- inla.mesh.2d(boundary=list(inla.mesh.boundary(mesh)),\n",
"                         max.edge=", min(mesh.edgelengths())/2,",
                         min.angle=21,
                         max.n=64000,
                         cutoff=0)\n"
                    )
            }
            out
        })
        mesh.proj <- reactive({
            if (is.null(mesh())) {
                NULL
            } else {
                INLA::inla.mesh.projector(mesh(), dims=c(500, 500))
            }
        })
        finemesh.proj <- reactive({
            if (is.null(finemesh())) {
                NULL
            } else {
                INLA::inla.mesh.projector(finemesh(), lattice=mesh.proj()$lattice)
            }
        })
        mesh.spde <- reactive({
            message("mesh.spde")
            if ((input$panel == "Input") || !input$assess || is.null(mesh())) {
                NULL
            } else {
            message("mesh.spde!")
                INLA::inla.spde2.pcmatern(mesh(),
                                          alpha=input$corr.nu+1,
                                          prior.range=c(1, 0.5),
                                          prior.sigma=c(1, 0.5))
            }
        })
        mesh.Q <- reactive({
            message("mesh.Q")
            if (!input$assess || is.null(mesh.spde())) {
                NULL
            } else {
            message("mesh.Q!")
                INLA::inla.spde.precision(mesh.spde(), theta=log(c(input$corr.range, 1)))
            }
        })
        mesh.S <- reactive({
            message("mesh.S")
            if (!input$assess || is.null(mesh.Q())) {
                NULL
            } else {
                message("mesh.S!")
                INLA::inla.qinv(mesh.Q(), reordering=INLA::inla.reorderings())
            }
        })
        finemesh.spde <- reactive({
            message("finemesh.spde")
            if ((input$panel == "Input") || !input$assess || is.null(finemesh())) {
                NULL
            } else {
                message("finemesh.spde!")
                INLA::inla.spde2.pcmatern(finemesh(),
                                          alpha=input$corr.nu+1,
                                          prior.range=c(1, 0.5),
                                          prior.sigma=c(1, 0.5))
            }
        })
        finemesh.Q <- reactive({
            message("finemesh.Q")
            if (!input$assess || is.null(finemesh.spde())) {
                NULL
            } else {
            message("finemesh.Q!")
                INLA::inla.spde.precision(finemesh.spde(), theta=log(c(input$corr.range, 1)))
            }
        })
        finemesh.S <- reactive({
            message("finemesh.S")
            if (!input$assess || is.null(finemesh.Q())) {
                NULL
            } else {
            message("finemesh.S!")
                INLA::inla.qinv(finemesh.Q(), reordering=INLA::inla.reorderings())
            }
        })
        observeEvent(c(input$overlay, input$sd.bound), {
            if (("corr" %in% input$overlay) && input$sd.bound) {
                updateCheckboxInput(session, "sd.bound", value=FALSE)
            }
        })
        mesh.sd <- reactive({
            message("mesh.sd")
            if (!input$assess || is.null(mesh.S())) {
                NULL
            } else {
                message("mesh.sd!")
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
            message("finemesh.sd")
            if (!input$assess || is.null(finemesh.S())) {
                NULL
            } else {
                message("finemesh.sd!")
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
        mesh.corr <- reactive({
            message("mesh.corr")
            message("mesh.corr A")
            A <- A.click()
            message("mesh.corr A,s")
            s <- SD.information()
            message("mesh.corr, s")
            message(paste("mesh.corr:", is.null(A), is.null(s)))
            if (!input$assess || !("corr" %in% input$overlay) ||
                is.null(mesh.sd()) || is.null(A) || is.null(s)) {
                NULL
            } else {
                message("mesh.corr!")
                proj <- mesh.proj()
                corr <- proj$proj$A %*% as.vector(inla.qsolve(mesh.Q(), t(A$mesh)))
                matrix(corr, length(proj$x), length(proj$y)) / mesh.sd() / s["SD","Mesh"]
            }
        })
        finemesh.corr <- reactive({
            message("finemesh.corr")
            A <- A.click()
            s <- SD.information()
            if (!input$assess || !("corr" %in% input$overlay) ||
                is.null(finemesh.sd()) || is.null(A) || is.null(s)) {
                NULL
            } else {
                message("finemesh.corr!")
                proj <- finemesh.proj()
                corr <- proj$proj$A %*% inla.qsolve(finemesh.Q(), t(A$finemesh))
                matrix(corr, length(proj$x), length(proj$y)) / finemesh.sd() / s["SD","Fine"]
            }
        })

        update.limits <- observe({
            lim <- 2^ceiling(log2(max(input$offset[2]*1.2,
                                      max(diff(input.xlim()), diff(input.ylim()))*1)))
            updateSliderInput(session, "offset",
                              min=lim/100, max=lim, step=lim/100)

            lim <- 2^ceiling(log2(max(input$max.edge[2]*1.2,
                                      max(diff(input.xlim()), diff(input.ylim()))*1)))
            updateSliderInput(session, "max.edge",
                              min=lim/100, max=lim, step=lim/100)

            lim <- 2^ceiling(log2(max(input$cutoff*1.2, input$max.edge[1]*1.2)))
            updateSliderInput(session, "cutoff",
                              min=lim/100, max=lim, step=lim/100)

            if (input$panel == "Display") {
                lim <- 2^ceiling(log2(max(input$corr.range*1.2,
                                          input$max.edge[1]*12,
                                          input$max.edge[2]*1.2)))
                updateSliderInput(session, "corr.range",
                                  min=lim/100, max=lim, step=lim/100)
            }
})
        
        meta <- reactive({
            if (input$panel != "Display") {
                m <- matrix(NA, 2, 3)
            } else {
                m1 <- mesh()
                if (is.null(m1)) {
                    m <- cbind(NA, NA, NA)
                } else {
                    m <- rbind(c(m1$n, range(mesh.edgelengths())))
                }
                m2 <- finemesh()
                if (is.null(m2)) {
                    m <- rbind(m, c(NA, NA, NA))
                } else {
                    m <- rbind(m, c(m2$n, range(finemesh.edgelengths())))
                }
            }
            rownames(m) <- c("Mesh", "Fine")
            colnames(m) <- c("nV", "minE", "maxE")
            m <- as.data.frame(m)
            m$nV <- as.integer(m$nV)
            m
        })
        
        meshplot.expr <- quote({
            message("meshplot.expr")
            if (input$panel != "Display") {
                return()
            }
            if (!input$assess && !is.null(mesh())) {
                message("meshplot.expr: basic")
                plot(NA, type="n", xlim=xlim(), ylim=ylim(),
                     asp=1, main="", xlab="Easting", ylab="Northing")
                plot(mesh(), add=TRUE)
            } else if (input$assess && !is.null(finemesh())) {
                message("meshplot.expr: assess")
                co <- mesh.corr()
                message(paste("meshplot.expr: mesh.corr", is.null(co)))
                
                if (("corr" %in% input$overlay) && !is.null(co)) {
                    message(paste("meshplot.expr: corr", input$assessment))
                    col <- colorRampPalette(c("blue", "white", "red"))
                    n.col <- 1+16
                    c.m <- co
                    c.f <- finemesh.corr()
                    c.r <- c.m - c.f
                    zlim.abs <- c(-1, 1)
                    zlim.rel <- c(-1, 1) * max(abs(c.r), na.rm=TRUE)
                    if (input$assessment == "coarse") {
                        message("meshplot.expr: coarse")
                        fields::image.plot(mesh.proj()$x, mesh.proj()$y, c.m,
                                           zlim=zlim.abs,
                                           col=col(n.col), asp=1, main="",
                                           xlab="Easting", ylab="Northing")
                    } else if (input$assessment == "rel") {
                        message("meshplot.expr: rel")
                        fields::image.plot(mesh.proj()$x, mesh.proj()$y, c.r,
                                           zlim=zlim.rel,
                                           col=col(n.col), asp=1, main="",
                                           xlab="Easting", ylab="Northing")
                    } else if (input$assessment == "fine") {
                        message("meshplot.expr: fine")
                        fields::image.plot(mesh.proj()$x, mesh.proj()$y, c.f,
                                           zlim=zlim.abs,
                                           col=col(n.col), asp=1, main="",
                                           xlab="Easting", ylab="Northing")
                    }
                } else {
                    message(paste("meshplot.expr: sd", input$assessment))
                    col <- colorRampPalette(c("blue", "white", "red"))
                    n.col <- 1+16
                    sd.m <- mesh.sd()
                    sd.f <- finemesh.sd()
                    sd.r <- sd.m / sd.f
                    zlim.abs <- range(c(range(sd.m, na.rm=TRUE),
                                        range(sd.f, na.rm=TRUE)))
                    zlim.abs <- exp(max(abs(log(zlim.abs))))
                    zlim.abs <- 1 + c(-1,1)*(zlim.abs-1)
                    zlim.rel <- range(sd.r, na.rm=TRUE)
                    zlim.rel <- exp(max(abs(log(zlim.rel))))
                    zlim.rel <- 1 + c(-1,1)*min(2, max(0.1, zlim.rel-1))
                    zlim.abs <- zlim.rel <- c(0.5, 1.5)
                    sd.m <- matrix(pmin(1.5, pmax(0.5, as.vector(sd.m))), nrow(sd.m), ncol(sd.m))
                    sd.f <- matrix(pmin(1.5, pmax(0.5, as.vector(sd.f))), nrow(sd.f), ncol(sd.f))
                    sd.r <- matrix(pmin(1.5, pmax(0.5, as.vector(sd.r))), nrow(sd.r), ncol(sd.r))
                    if (input$assessment == "coarse") {
                        message("meshplot.expr: coarse")
                        fields::image.plot(mesh.proj()$x, mesh.proj()$y, sd.m,
                                           zlim=zlim.abs,
                                           col=col(n.col), asp=1, main="",
                                           xlab="Easting", ylab="Northing")
                    } else if (input$assessment == "rel") {
                        message("meshplot.expr: rel")
                        fields::image.plot(mesh.proj()$x, mesh.proj()$y, sd.r,
                                           zlim=zlim.rel,
                                           col=col(n.col), asp=1, main="",
                                           xlab="Easting", ylab="Northing")
                    } else if (input$assessment == "fine") {
                        message("meshplot.expr: fine")
                        fields::image.plot(mesh.proj()$x, mesh.proj()$y, sd.f,
                                           zlim=zlim.abs,
                                           col=col(n.col), asp=1, main="",
                                           xlab="Easting", ylab="Northing")
                    }
                }
                if ("mesh" %in% input$overlay) {
                    message("meshplot.expr: mesh overlay")
                    plot(mesh(), add=TRUE)
                }
                if ("finemesh" %in% input$overlay) {
                    message("meshplot.expr: fine overlay")
                    plot(finemesh(), add=TRUE)
                }
            }
        })

        output$meshplot <-
            renderPlot(meshplot.expr, quoted=TRUE, height=1000, units="px")

        userinput.xlim <- reactive({
            lim1 <- if (is.null(boundary.loc.input())) NA else range(coordinates(boundary.loc.input())[,1], na.rm=TRUE)
            lim2 <- if (is.null(mesh.loc.input())) NA else range(coordinates(mesh.loc.input())[,1], na.rm=TRUE)
            r <- c(lim1, lim2)
            if (all(is.na(r))) {
                r <- c(0, 1)
            } else {
                r <- range(r, na.rm=TRUE)
            }
            r
        })
        userinput.ylim <- reactive({
            lim1 <- if (is.null(boundary.loc.input())) NA else range(coordinates(boundary.loc.input())[,2], na.rm=TRUE)
            lim2 <- if (is.null(mesh.loc.input())) NA else range(coordinates(mesh.loc.input())[,2], na.rm=TRUE)
            r <- c(lim1, lim2)
            if (all(is.na(r))) {
                r <- c(0, 1)
            } else {
                r <- range(r, na.rm=TRUE)
            }
            r
        })

        input.xlim <- reactive({
            lim1 <- userinput.xlim()
            lim2 <- if (length(input$loc.usage) > 0) {
                if (is.null(loc())) NA else range(coordinates(loc())[,1], na.rm=TRUE)
            } else {
                NA
            }
            r <- c(lim1, lim2)
            if (all(is.na(r))) {
                r <- c(0, 1)
            } else {
                r <- range(r, na.rm=TRUE)
            }
            r
        })
        input.ylim <- reactive({
            lim1 <- userinput.ylim()
            lim2 <- if (length(input$loc.usage) > 0) {
                if (is.null(loc())) NA else range(coordinates(loc())[,2], na.rm=TRUE)
            } else {
                NA
            }
            r <- c(lim1, lim2)
            if (all(is.na(r))) {
                r <- c(0, 1)
            } else {
                r <- range(r, na.rm=TRUE)
            }
            r
        })

        xlim <- reactive({
            lim1 <- if (is.null(mesh())) NA else range(mesh()$loc[,1], na.rm=TRUE)
            r <- c(lim1)
            if (all(is.na(r))) {
                r <- c(0, 1)
            } else {
                r <- range(r, na.rm=TRUE)
            }
            r
        })
        ylim <- reactive({
            lim1 <- if (is.null(mesh())) NA else range(mesh()$loc[,2], na.rm=TRUE)
            r <- lim1
            if (all(is.na(r))) {
                r <- c(0, 1)
            } else {
                r <- range(r, na.rm=TRUE)
            }
            r
        })

        inputplot.expr <- quote({
            if (input$panel != "Input") {
                return()
            }
            plot(NA, type="n", main="", xlim=input.xlim(), ylim=input.ylim(),
                 xlab="Easting", ylab="Northing")
            points(mesh.loc(), col=2)
            points(boundary.loc(), col=4)
        })

        output$inputplot <-
            renderPlot(inputplot.expr, quoted=TRUE, height=800, units="px")

        output$meta <- renderTable({meta()}, digits=3, rownames=TRUE)

        clicks <- reactiveValues(loc.click=NULL)
        observeEvent(input$meshplot.click, {
            message("loc.click")
            if (is.null(input$meshplot.click)) {
                clicks$loc.click <- NULL
            } else {
                clicks$loc.click <- cbind(input$meshplot.click$x, input$meshplot.click$y)
            }
            message(paste("loc.click: loc.click =", paste(clicks$loc.click, collapse=", ")))
        })
        A.click <- reactive({
            message(paste("A.click:",
                          paste0("loc.click = ", paste(clicks$loc.click, collapse=", ")),
                          (input$panel == "Input"), !input$assess,
                          is.null(mesh()), is.null(finemesh()), (is.null(clicks$loc.click))))
            if ((input$panel == "Input") || !input$assess ||
                is.null(mesh()) || is.null(finemesh()) || (is.null(clicks$loc.click))) {
                NULL
            } else {
                message("A.click!")
                A.mesh <- INLA::inla.spde.make.A(mesh(), clicks$loc.click)
                if (!(sum(A.mesh) > 0)) {
                    return(NULL)
                }
                A.finemesh <- INLA::inla.spde.make.A(finemesh(), clicks$loc.click)
                if (!(sum(A.finemesh) > 0)) {
                    return(NULL)
                }
                message("A.click!!")
                list(mesh=A.mesh, finemesh=A.finemesh)
            }
        })
        SD.information <- reactive({
            message(paste("SD.information",
                          !input$assess, is.null(mesh()), is.null(finemesh()),
                          is.null(A.click())))
            A <- A.click()
            if (!input$assess || is.null(mesh()) || is.null(finemesh()) || is.null(A)) {
                NULL
            } else {
                message("SD.information!")
                if (is.null(A)) {
                    out <- rbind(cbind(NA, NA, NA), cbind(NA, NA, NA))
                } else {
                    message("SD.information!!")
                    sd.bound.mesh <- as.vector(A$mesh %*% diag(mesh.S())^0.5)
                    sd.bound.finemesh <- as.vector(A$finemesh %*% diag(finemesh.S())^0.5)
                    sd.mesh <- Matrix::rowSums(A$mesh * (A$mesh %*% mesh.S()))^0.5
                    sd.finemesh <- Matrix::rowSums(A$finemesh * (A$finemesh %*% finemesh.S()))^0.5
                    out <- rbind(cbind(sd.mesh, sd.mesh/sd.finemesh, sd.finemesh),
                                 cbind(sd.bound.mesh,
                                       sd.bound.mesh/sd.bound.finemesh,
                                       sd.bound.finemesh))
                }
                rownames(out) <- c("SD", "SD bound")
                colnames(out) <- c("Mesh", "Ratio", "Fine")
                out <- as.data.frame(out)
                out
            }
        })
        corr.information <- reactive({
            message(paste("corr.information",
                          !input$assess, !("corr" %in% input$overlay), is.null(mesh.corr()),
                          is.null(finemesh.corr())))
            if (!input$assess || !("corr" %in% input$overlay) || is.null(mesh.corr()) ||
                is.null(finemesh.corr())) {
                NULL
            } else {
                message("corr.information!")
                val <- as.vector(mesh.corr() - finemesh.corr())
                out <- matrix(NA, 1,5)
                out[1,1] <- min(val, na.rm=TRUE)
                out[1,2:4] <- quantile(val, c(0.25, 0.5, 0.75), na.rm=TRUE)
                out[1,5] <- max(val, na.rm=TRUE)
                rownames(out) <- c("Correlation error")
                colnames(out) <- c("min", "25%", "50%", "75%", "max")
                out <- as.data.frame(out)
                out
            }
        })
        
        output$SD.information <- renderTable({SD.information()},
                                             digits=3, rownames=TRUE)
        output$corr.information <- renderTable({corr.information()},
                                               digits=3, rownames=TRUE)

        output$code <- renderText({
            out <- "## Generated by meshbuilder()\n\n"
            spacing <- ""
            if (length(input$loc.usage) > 0) {
                out <- paste0(out,
                              "## Prepare seed points:\n",
                              attr(loc(), "code"))
                spacing <- "\n"
            }
            if (!is.null(boundary())) {
                out <- paste0(out, spacing,
                              "## Build boundary information:\n",
                              attr(boundary.loc(), "code"),
                              attr(boundary(), "code"))
                spacing <- "\n"
            }
            out <- paste0(out, spacing,
                          "## Build the mesh:\n",
                          attr(mesh.loc(), "code"),
                          attr(mesh(), "code"))
            spacing <- "\n"
            out <- paste0(out, spacing,
                          "## Plot the mesh:\n",
                          "plot(mesh)\n"
                          )
            if (input$assess) {
                out <- paste0(out, spacing, spacing,
                              "## Build a fine scale mesh for comparisons:\n",
                              attr(finemesh(), "code"), "\n")
            }
            out
        })
        
        output$crs.strings <- renderText({
            crs <- c()
            sp <- mesh.loc()
            if (!is.null(sp) && !is.na(proj4string(sp))) {
                crs <- c(crs, proj4string(sp))
            }
            sp <- boundary.loc()
            if (!is.null(sp) && !is.na(proj4string(sp))) {
                crs <- c(crs, proj4string(sp))
            }
            crs <- unique(crs)
            out <- paste0(crs, collapse="\n")
            out
        })

        output$input.coordinates <- renderText({
            xy_str <- function(e) {
                if (is.null(e)) {
                    "NULL\n"
                } else {
                    paste0("x=", e$x, ", y=", e$y, "\n")
                }
            }
            paste0(xy_str(input$inputplot.hover))
        })
        output$loc.click.coordinates <- renderText({
            xy_str <- function(e) {
                if (is.null(e)) {
                    "NULL\n"
                } else {
                    paste0("x=", e[1,1], ", y=", e[1,2], "\n")
                }
            }
            paste0(xy_str(clicks$loc.click))
        })

        output$debug <- renderText({
            out <- paste0(input$boundary.loc.name)
            out <- paste0(out, "\n", input$mesh.loc.name)
            out
        })
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
