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

spatial.object.choices <- function() {
    valid.classes <- c("matrix",
                       "inla.mesh.segment",
                       "SpatialPoints", "SpatialPointsDataFrame",
                       "SpatialLines", "SpatialLinesDataFrame",
                       "SpatialPolygons", "SpatialPolygonsDataFrame")
    obj <- with(globalenv(), ls())
    ok <- vapply(obj,
                 function(name) {
                     cl <- class(globalenv()[[name]])
                     (length(cl) > 0) && (length(intersect(cl, valid.classes)) > 0)
                 },
                 TRUE)
    obj[ok]
}

pretty.axis.info <- function(lim, value=NA) {
    message(paste("pretty.axis.info <-",
                  "lim =", paste(lim, collapse=", "),
                  "val = ", paste(value, collapse=", ")))
    limits <- range(c(lim, value), na.rm=TRUE)
    maxi <- c(1.5, 2, 3, 4, 6, 8, 10)
    cutoff <- maxi * 9/10
    step <- c(1e-2, 1e-2, 1e-2, 2e-2, 2e-2, 5e-2, 5e-2)
    mini <- step
    maxi <- c(maxi, maxi*10)
    step <- c(step, step*10)
    mini <- c(mini, mini*10)
    cutoff <- c(cutoff, cutoff*10)
    level <- floor(log10(limits[2]))
    maxi <- maxi * 10^level
    mini <- mini * 10^level
    step <- step * 10^level
    cutoff <- cutoff * 10^level
    ## Find smallest k such that limits[2] < maxi[k]
    k <- min(which(limits[2] < maxi))
    ## Move to next k if limits[2] > cutoff[k]
    values <- mini[k] + round((value - mini[k]) / step[k]) * step[k]
    limits <- range(c(lim, value, values), na.rm=TRUE)
    while (limits[2] > cutoff[k]) {
        k <- k+1
        values <- mini[k] + round((value - mini[k]) / step[k]) * step[k]
        limits <- range(c(lim, value, value), na.rm=TRUE)
    }

    message(paste("pretty.axis.info ->",
                  "lim =", paste(c(mini[k], maxi[k]), collapse=", ")))
    list(lim=c(mini[k], maxi[k]), step=step[k])
}

pretty.axis.value <- function(axis.info, value) {
    pmax(axis.info$lim[1], pmin(axis.info$lim[2],
                                axis.info$lim[1] + round((value - axis.info$lim[1]) /
                                                         axis.info$step) * axis.info$step))
}

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
                           uiOutput("overlay.ui"),
                           uiOutput("assess.resolution.ui")
                           ),
                    column(6,
                           checkboxInput("assess.advanced", label = "Advanced", value=FALSE),
                           uiOutput("assess.quantity.ui")
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
                                                    choices=spatial.object.choices()),
                                        selectInput("boundary.loc.name.usage",
                                                    label=NULL,
                                                    multiple=TRUE,
                                                    choices=c("points", "domain")),
                                        selectInput("mesh.loc.name",
                                                    label="Mesh vertex seed points variable name(s)",
                                                    multiple=TRUE,
                                                    choices=spatial.object.choices()),
                                        h4("User defined boundaries"),
                                        selectInput("boundary1.name",
                                                    label="Inner boundary variable name(s)",
                                                    multiple=TRUE,
                                                    choices=with(globalenv(), ls())),
                                        selectInput("boundary2.name",
                                                    label="Outer boundary variable name(s)",
                                                    multiple=TRUE,
                                                    choices=with(globalenv(), ls())),
                                        hr(),
                                        textInput("crs.mesh", label="Mesh CRS"),
                                        verbatimTextOutput("crs.strings")
                                        ),
                                 column(8,
                                        plotOutput("inputplot", height="auto",
                                                   hover="inputplot.hover")
                                        )
                             ),
                             hr(),
                             verbatimTextOutput("input.loc.coordinates")
                             ),
                    tabPanel("Display",
                             plotOutput("meshplot", height="auto",
                                        click="meshplot.click"
                                        ),
                             fluidRow(
                                 column(3,
                                        tableOutput("meshmeta")
                                        ),
                                 column(3,
                                        tableOutput("click.information")
                                        ),
                                 column(6,
                                        tableOutput("field.information")
                                        )
                             ),
                             verbatimTextOutput("meshplot.loc.coordinates")
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
                metadata$mesh <- NULL
            } else {
                try(message(paste("crs.mesh               =",
                                  inla.CRSargs(inla.CRS(input$crs.mesh)))))
                try(message(paste("CRS(mesh.loc)          =", proj4string(loc1))))
                try(message(paste("CRS(boundary[[1]]$crs) =", inla.CRSargs(bnd1[[1]]$crs))))
                try(message(paste("CRS(boundary[[2]]$crs) =", inla.CRSargs(bnd1[[2]]$crs))))
                out <- INLA::inla.mesh.2d(loc=loc1,
                                          boundary=bnd1,
                                          max.edge=input$max.edge,
                                          min.angle=rev(input$min.angle),
                                          max.n=c(48000, 16000),
                                          cutoff=input$cutoff,
                                          crs=inla.CRS(input$crs.mesh))
                attr(out, "code") <- paste0(
"mesh <- inla.mesh.2d(", loc.code, bnd.code,
                    "max.edge=c(", input$max.edge[1], ", ", input$max.edge[2], "),
                     min.angle=c(", input$min.angle[2], ", ", input$min.angle[1], "),
                     max.n=c(48000, 16000)
                     cutoff=", input$cutoff, ")\n")
                el <- sqrt(c(Matrix::rowSums((out$loc[out$graph$tv[,2],] -
                                              out$loc[out$graph$tv[,1],])^2),
                             Matrix::rowSums((out$loc[out$graph$tv[,3],] -
                                              out$loc[out$graph$tv[,2],])^2),
                             Matrix::rowSums((out$loc[out$graph$tv[,1],] -
                                              out$loc[out$graph$tv[,3],])^2)))
                metadata$mesh <- cbind(out$n, min(el), max(el))
            }
            metadata$fine <- NULL
            out
        })
        fine <- reactive({
            if (is.null(mesh())) {
                out <- NULL
                metadata$fine <- NULL
            } else {
                message("fine")
                out <- INLA::inla.mesh.2d(loc=mesh()$loc,
                                          boundary=list(INLA::inla.mesh.boundary(mesh())),
                                          max.edge=min(mesh.edgelengths())/2,
                                          min.angle=21,
                                          max.n=64000,
                                          cutoff=0)
                attr(out, "code") <- paste0(
"fine <- inla.mesh.2d(boundary=list(inla.mesh.boundary(mesh)),\n",
"                         max.edge=", min(mesh.edgelengths())/2,",
                         min.angle=21,
                         max.n=64000,
                         cutoff=0)\n"
                    )
                el <- sqrt(c(Matrix::rowSums((out$loc[out$graph$tv[,2],] -
                                              out$loc[out$graph$tv[,1],])^2),
                             Matrix::rowSums((out$loc[out$graph$tv[,3],] -
                                              out$loc[out$graph$tv[,2],])^2),
                             Matrix::rowSums((out$loc[out$graph$tv[,1],] -
                                              out$loc[out$graph$tv[,3],])^2)))
                metadata$fine <- cbind(out$n, min(el), max(el))
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
        fine.edgelengths <- reactive({
            m <- fine()
            if (is.null(m)) {
                c(NA)
            } else {
                sqrt(c(Matrix::rowSums((m$loc[m$graph$tv[,2],] - m$loc[m$graph$tv[,1],])^2),
                       Matrix::rowSums((m$loc[m$graph$tv[,3],] - m$loc[m$graph$tv[,2],])^2),
                       Matrix::rowSums((m$loc[m$graph$tv[,1],] - m$loc[m$graph$tv[,3],])^2)))
            }
        })

        mesh.proj <- reactive({
            if (is.null(mesh())) {
                NULL
            } else {
                INLA::inla.mesh.projector(mesh(), dims=c(500, 500))
            }
        })
        fine.proj <- reactive({
            if (is.null(fine())) {
                NULL
            } else {
                INLA::inla.mesh.projector(fine(), lattice=mesh.proj()$lattice)
            }
        })
        mesh.spde <- reactive({
            message("mesh.spde")
            if (is.null(mesh())) {
                NULL
            } else {
                message("mesh.spde!")
                INLA::inla.spde2.pcmatern(mesh(),
                                          alpha=input$corr.nu+1,
                                          prior.range=c(1, 0.5),
                                          prior.sigma=c(1, 0.5))
            }
        })
        fine.spde <- reactive({
            message("fine.spde")
            if (is.null(mesh())) {
                NULL
            } else {
                message("fine.spde!")
                INLA::inla.spde2.pcmatern(fine(),
                                          alpha=input$corr.nu+1,
                                          prior.range=c(1, 0.5),
                                          prior.sigma=c(1, 0.5))
            }
        })
        mesh.Q <- reactive({
            message("mesh.Q")
            if (is.null(mesh.spde())) {
                NULL
            } else {
            message("mesh.Q!")
                INLA::inla.spde.precision(mesh.spde(), theta=log(c(input$corr.range, 1)))
            }
        })
        fine.Q <- reactive({
            message("fine.Q")
            if (is.null(fine.spde())) {
                NULL
            } else {
            message("fine.Q!")
                INLA::inla.spde.precision(fine.spde(), theta=log(c(input$corr.range, 1)))
            }
        })
        mesh.S <- reactive({
            message("mesh.S")
            if (is.null(mesh.Q())) {
                NULL
            } else {
                message("mesh.S!")
                INLA::inla.qinv(mesh.Q(), reordering=INLA::inla.reorderings())
            }
        })
        fine.S <- reactive({
            message("fine.S")
            if (is.null(fine.Q())) {
                NULL
            } else {
                message("fine.S!")
                INLA::inla.qinv(fine.Q(), reordering=INLA::inla.reorderings())
            }
        })
        mesh.sd <- reactive({
            message("mesh.sd")
            if (is.null(mesh.S())) {
                message("mesh.sd: NULL")
                NULL
            } else {
                message("mesh.sd!")
                proj <- mesh.proj()
                v <- Matrix::rowSums(proj$proj$A * (proj$proj$A %*% mesh.S()))
                v[!proj$proj$ok] <- NA
                matrix(v^0.5, length(proj$x), length(proj$y))
            }
        })
        fine.sd <- reactive({
            message("fine.sd")
            if (is.null(fine.S())) {
                NULL
            } else {
                message("fine.sd!")
                proj <- fine.proj()
                v <- Matrix::rowSums(proj$proj$A * (proj$proj$A %*% fine.S()))
                v[!proj$proj$ok] <- NA
                matrix(v^0.5, length(proj$x), length(proj$y))
            }
        })
        mesh.sd.deviation.approx <- reactive({
            message("mesh.sd.deviation.approx")
            if (is.null(mesh.sd()) || is.null(mesh.proj())) {
                message("mesh.sd.deviation.approx: NULL")
                NULL
            } else {
                val0 <- mesh.sd()
                proj <- mesh.proj()
                val <- proj$proj$A %*% (
                    as.vector(t(proj$proj$A[proj$proj$ok,,drop=FALSE]) %*%
                              as.vector(val0)[proj$proj$ok]) /
                    colSums(proj$proj$A[proj$proj$ok,,drop=FALSE]))
                val[!proj$proj$ok] <- NA
                matrix(1 + (as.vector(val0)-val),
                       length(proj$x), length(proj$y))
            }
        })
        fine.sd.deviation.approx <- reactive({
            message("fine.sd.deviation.approx")
            if (is.null(fine.sd()) || is.null(fine.proj())) {
                message("fine.sd.deviation.approx: NULL")
                NULL
            } else {
                val0 <- fine.sd()
                proj <- fine.proj()
                val <- proj$proj$A %*% (
                    as.vector(t(proj$proj$A[proj$proj$ok,,drop=FALSE]) %*%
                              as.vector(val0)[proj$proj$ok]) /
                    colSums(proj$proj$A[proj$proj$ok,,drop=FALSE]))
                val[!proj$proj$ok] <- NA
                matrix(1 + (as.vector(val0)-val),
                       length(proj$x), length(proj$y))
            }
        })

        mesh.sd.bound <- reactive({
            message("mesh.sd.bound")
            if (is.null(mesh.S())) {
                NULL
            } else {
                message("mesh.sd.bound!")
                proj <- mesh.proj()
                INLA::inla.mesh.project(proj, field=diag(mesh.S())^0.5)
            }
        })
        fine.sd.bound <- reactive({
            message("fine.sd.bound")
            if (is.null(fine.S())) {
                NULL
            } else {
                message("fine.sd.bound!")
                proj <- fine.proj()
                INLA::inla.mesh.project(proj, field=diag(fine.S())^0.5)
            }
        })
        mesh.corr <- reactive({
            message("mesh.corr")
            if (is.null(mesh.sd()) || is.null(meshplot.A.mesh()) ||
                is.null(click.information()) || is.null(mesh.proj())) {
                message("mesh.corr: NULL")
                NULL
            } else {
                message("mesh.corr!")
                A <- meshplot.A.mesh()
                s <- click.information()
                proj <- mesh.proj()
                print(str(A))
                print(str(s))
                print(str(proj))
                corr <- proj$proj$A %*% as.vector(inla.qsolve(mesh.Q(), t(A)))
                out <- matrix(corr, length(proj$x), length(proj$y)) / mesh.sd() / s["SD","Mesh"]
                message("mesh.corr: end")
                out
            }
        })
        fine.corr <- reactive({
            message("fine.corr")
            if (is.null(fine.sd()) || is.null(meshplot.A.fine()) ||
                is.null(click.information) || is.null(fine.proj())) {
                NULL
            } else {
                message("fine.corr!")
                A <- meshplot.A.fine()
                s <- click.information()
                proj <- fine.proj()
                corr <- proj$proj$A %*% inla.qsolve(fine.Q(), t(A))
                matrix(corr, length(proj$x), length(proj$y)) / fine.sd() / s["SD","Fine"]
            }
        })

        axis.infos <- reactiveValues(offset=NULL, max.edge=NULL, cutoff=NULL, corr.range=NULL)
        axis.update <- reactiveValues(offset=FALSE, max.edge=FALSE, cutoff=FALSE, corr.range=NULL)
        observeEvent(c(input$offset, limits$input.xlim, limits$input$ylim), {
            message("Observe offset limit recalculation")
            val <- input$offset
            new.info <- pretty.axis.info(range(c(val,
                                                 diff(limits$input.xlim)/5,
                                                 diff(limits$input.ylim)/5)),
                                         val)
            print(axis.infos$offset$lim)
            print(new.info$lim)
            axis.update$offset <- !identical(axis.infos$offset, new.info)
            message(paste("Update offset", axis.update$offset))
            if (axis.update$offset) {
                axis.infos$offset <- new.info
            }
        })
        observeEvent(c(input$max.edge, limits$input.xlim, limits$input.ylim), {
            message("Observe max.edge limit recalculation")
            val <- input$max.edge
            new.info <- pretty.axis.info(range(c(val,
                                                 diff(limits$input.xlim)/5,
                                                 diff(limits$input.ylim)/5)),
                                         val)
            axis.update$max.edge <- !identical(axis.infos$max.edge, new.info)
            if (axis.update$max.edge) {
                axis.infos$max.edge <- new.info
            }
        })
        observeEvent(c(input$cutoff, input$max.edge), {
            message("Observe cutoff limit recalculation")
            val <- min(input$cutoff, input$max.edge[1])
            new.info <- pretty.axis.info(range(c(val, input$max.edge[1])),
                                         val)
            axis.update$cutoff <- !identical(axis.infos$cutoff, new.info)
            if (axis.update$cutoff) {
                axis.infos$cutoff <- new.info
            }
        })
        observeEvent(c(input$corr.range, limits$input.xlim, limits$input.ylim), {
            message("Observe corr.range limit recalculation")
            val <- input$corr.range
            new.info <- pretty.axis.info(range(c(val,
                                                 diff(limits$input.xlim),
                                                 diff(limits$input.ylim))),
                                         val)
            axis.update$corr.range <- !identical(axis.infos$corr.range, new.info)
            if (axis.update$corr.range) {
                axis.infos$corr.range <- new.info
            }
        })
        observeEvent(axis.update$offset, {
            message(paste("Observe offset limit update", axis.update$offset))
            if (axis.update$offset) {
                axis.update$offset <- FALSE
                info <- axis.infos$offset
                new.val <- pretty.axis.value(info, input$offset)
                print(paste("Old values", paste(input$offset, collapse=", ")))
                print(paste("New limits", paste(info$lim, collapse=", ")))
                print(paste("New values", paste(new.val, collapse=", ")))
                updateSliderInput(session, "offset",
                                  min=info$lim[1], max=info$lim[2],
                                  step=info$step, value=new.val)
            }
        }, priority=10)
        observeEvent(axis.update$max.edge, {
            message(paste("Observe max.edge limit update", axis.update$max.edge))
            if (axis.update$max.edge) {
                axis.update$max.edge <- FALSE
                info <- axis.infos$max.edge
                new.val <- pretty.axis.value(info, input$max.edge)
                print(paste("Old values", paste(input$max.edge, collapse=", ")))
                print(paste("New limits", paste(info$lim, collapse=", ")))
                print(paste("New values", paste(new.val, collapse=", ")))
                updateSliderInput(session, "max.edge",
                                  min=info$lim[1], max=info$lim[2],
                                  step=info$step, value=new.val)
            }
        }, priority=9)
        observeEvent(c(axis.update$cutoff, input$cutoff, input$max.edge), {
            message(paste("Observe cutoff limit update", axis.update$cutoff))
            if (axis.update$cutoff || (input$cutoff > input$max.edge[1])) {
                axis.update$cutoff <- FALSE
                info <- axis.infos$cutoff
                new.val <- pretty.axis.value(info, min(input$cutoff, input$max.edge[1]))
                print(paste("Old values", paste(input$cutoff, collapse=", ")))
                print(paste("New limits", paste(info$lim, collapse=", ")))
                print(paste("New values", paste(new.val, collapse=", ")))
                updateSliderInput(session, "cutoff",
                                  min=info$lim[1], max=info$lim[2],
                                  step=info$step, value=new.val)
            }
        }, priority=8)
        observeEvent(axis.update$corr.range, {
            message(paste("Observe corr.range limit update", axis.update$corr.range))
            if (axis.update$corr.range) {
                axis.update$corr.range <- FALSE
                info <- axis.infos$corr.range
                new.val <- pretty.axis.value(info, input$corr.range)
                print(paste("Old values", paste(input$corr.range, collapse=", ")))
                print(paste("New limits", paste(info$lim, collapse=", ")))
                print(paste("New values", paste(new.val, collapse=", ")))
                updateSliderInput(session, "corr.range",
                                  min=info$lim[1], max=info$lim[2],
                                  step=info$step, value=new.val)
            }
        }, priority=7)
        
        metadata <- reactiveValues(mesh=NULL,
                                   fine=NULL)
        meshmeta <- reactive({
            m1 <- m2 <- NULL
            if (!is.null(metadata$mesh)) {
                m1 <- rbind(metadata$mesh)
                rownames(m1) <- "Mesh"
            }
            if (!is.null(metadata$fine)) {
                m2 <- rbind(metadata$fine)
                rownames(m2) <- "Fine"
            }
            out <- rbind(m1, m2)
            colnames(out) <- c("nV", "minE", "maxE")
            out <- as.data.frame(out)
            out$nV <- as.integer(out$nV)
            out
        })

        observeEvent(c(input$assess.quantity, input$assess.resolution), {
            if (!is.null(input$assess.quantity)) {
                if (input$assess.quantity == "deviation.approx") {
                    if (input$assess.resolution == "rel") {
                        updateRadioButtons(session, "assess.resolution", selected="mesh")
                    }
                } else
                if (input$assess.quantity == "deviation") {
                    updateRadioButtons(session, "assess.quantity", selected="sd")
                    updateRadioButtons(session, "assess.resolution", selected="rel")
                }
            }
        }, priority=10)

        output$assess.quantity.ui <- renderUI({
            sel <- isolate(input$assess.quantity)
            default <- 1
            choices <- list("SD" = "sd",
                            "Correlation" = "corr",
                            "Approximate SD deviation" = "deviation.approx",
                            "SD bound" = "sd.bound",
                            "SD/bound" = "sd/bound")
            if (input$assess.advanced) {
                choices <- c(choices, list("SD deviation" = "deviation"))
            }
            if (is.null(sel) || !(sel %in% choices)) {
                sel <- choices[[default]]
            }
            radioButtons("assess.quantity", label = "Quantity",
                         choices = choices, selected = sel)
        })
        output$overlay.ui <- renderUI({
            if (is.null(input$assess.advanced)) {
                return(NULL)
            }
            sel <- isolate(input$overlay)
            default <- 2
            if (input$assess.advanced) {
                choices <- list("None" = "none",
                                "Mesh" = "mesh",
                                "Fine" = "fine")
            } else {
                choices <- list("None" = "none",
                                "Mesh" = "mesh")
            }
            if (is.null(sel) || !(sel %in% choices)) {
                sel <- choices[[default]]
            }
            radioButtons("overlay", label="Overlay",
                         choices = choices, selected = sel)
        })
        output$assess.resolution.ui <- renderUI({
            if (is.null(input$assess.quantity)) {
                return(NULL)
            }
            sel <- isolate(input$assess.resolution)
            default <- 1
            if (input$assess.advanced) {
                if (input$assess.quantity == "corr") {
                    choices <- list("Mesh" = "mesh",
                                    "Fine" = "fine",
                                    "Mesh-Fine" = "rel")
                } else {
                    choices <- list("Mesh" = "mesh",
                                    "Fine" = "fine",
                                    "Mesh/Fine" = "rel")
                }
            } else {
                choices <- list("Mesh" = "mesh")
            }
            if (is.null(sel) || !(sel %in% choices)) {
                sel <- choices[[default]]
            }
            radioButtons("assess.resolution", label = "Resolution",
                         choices = choices, selected=sel)
            })
        
        meshplot.expr <- quote({
            message("meshplot.expr")
            col <- colorRampPalette(c("blue", "white", "red"))
            n.col <- 1+64
            if (!input$assess && !is.null(mesh())) {
                message("meshplot.expr: basic")
                zlim <- c(0.5, 1.5)
##                fields::image.plot(mesh.proj()$x, mesh.proj()$y,
##                                   matrix(0, length(mesh.proj()$x), length(mesh.proj()$y)),
##                                   zlim=zlim, 
##                                   xlim=xlim(), ylim=ylim(),
##                                   col=col(n.col), asp=1, main="",
##                                   xlab="Easting", ylab="Northing")
                fields::image.plot(range(mesh.proj()$x), range(mesh.proj()$y),
                                   matrix(0, 2, 2),
                                   zlim=zlim, 
                                   xlim=xlim(), ylim=ylim(),
                                   col=col(n.col), asp=1, main="",
                                   xlab="Easting", ylab="Northing")
                plot(mesh(), add=TRUE)
            } else if (input$assess) {
                message("meshplot.expr: assess")
                message(paste("meshplot.expr: assess.quantity", input$assess.quantity))
                message(paste("meshplot.expr: assess.resolution", input$assess.resolution))
                message(paste("meshplot.expr: overlay", input$overlay))
                zlim.sd <- c(0.5, 1.5)
                zlim.sd.rel <- c(0.5, 1.5)
                zlim.corr <- c(-1, 1)
                zlim.corr.rel <- c(-1, 1)/5
                zlim.sdbound <- c(1/sqrt(3), 1+(1-1/sqrt(3)))
                zlim.sdbound.rel <- c(0.5, 1.5)
                ## sd approximate deviation
                if (input$assess.quantity %in% c("deviation.approx")) {
                    zlim <- zlim.sd.rel
                    if (input$assess.resolution %in% c("mesh")) {
                        val <- mesh.sd.deviation.approx()
                    } else if (input$assess.resolution %in% c("fine")) {
                        val <- fine.sd.deviation.approx()
                    } else if (input$assess.resolution %in% c("rel")) {
                        val <- mesh.sd.deviation.approx() / fine.sd.deviation.approx()
                    }
                } else
                ## sd
                if (input$assess.quantity %in% c("sd")) {
                    zlim <- zlim.sd
                    if (input$assess.resolution %in% c("mesh")) {
                        val <- mesh.sd()
                    } else if (input$assess.resolution %in% c("fine")) {
                        val <- fine.sd()
                    } else if (input$assess.resolution %in% c("rel")) {
                        val <- mesh.sd() / fine.sd()
                        zlim <- zlim.sd.rel
                    }
                } else
                ## corr
                if (input$assess.quantity %in% c("corr")) {
                    zlim <- zlim.corr
                    if (input$assess.resolution %in% c("mesh")) {
                        val <- mesh.corr()
                    } else if (input$assess.resolution %in% c("fine")) {
                        val <- fine.corr()
                    } else if (input$assess.resolution %in% c("rel")) {
                        val <- mesh.corr()
                        if (!is.null(val)) {
                            val <- val - fine.corr()
                        }
                        zlim <- zlim.corr.rel
                    }
                } else
                ## sd.bound
                if (input$assess.quantity %in% c("sd.bound")) {
                    zlim <- zlim.sd
                    if (input$assess.resolution %in% c("mesh")) {
                        val <- mesh.sd.bound()
                    } else if (input$assess.resolution %in% c("fine")) {
                        val <- fine.sd.bound()
                    } else if (input$assess.resolution %in% c("rel")) {
                        val <- mesh.sd.bound() / fine.sd.bound()
                        zlim <- zlim.sd.rel
                    }
                } else
                ## sd/bound
                if (input$assess.quantity %in% c("sd/bound")) {
                    zlim <- zlim.sdbound
                    if (input$assess.resolution %in% c("mesh")) {
                        val <- mesh.sd() / mesh.sd.bound()
                    } else if (input$assess.resolution %in% c("fine")) {
                        val <- fine.sd() / fine.sd.bound()
                    } else if (input$assess.resolution %in% c("rel")) {
                        val <- (mesh.sd() / mesh.sd.bound()) / (fine.sd() / fine.sd.bound())
                        zlim <- zlim.sdbound.rel
                    }
                } else {
                    val <- NULL
                    zlim <- c(-1, 1)
                }
                if (is.null(val)) {
                    fields::image.plot(mesh.proj()$x, mesh.proj()$y,
                                       matrix(0, length(mesh.proj()$x), length(mesh.proj()$y)),
                                       zlim=zlim,
                                       xlim=xlim(), ylim=ylim(),
                                       col=col(n.col), asp=1, main="",
                                       xlab="Easting", ylab="Northing")
                } else {
                    val <- matrix(pmax(zlim[1], pmin(zlim[2], as.vector(val))),
                                  nrow(val), ncol(val))
                    fields::image.plot(mesh.proj()$x, mesh.proj()$y, val,
                                       zlim=zlim,
                                       xlim=xlim(), ylim=ylim(),
                                       col=col(n.col), asp=1, main="",
                                       xlab="Easting", ylab="Northing")
                }
                if (!is.null(clicks$meshplot.loc)) {
                    points(clicks$meshplot.loc, col=2)
                }
                if ("mesh" %in% input$overlay) {
                    message("meshplot.expr: mesh overlay")
                    plot(mesh(), add=TRUE)
                }
                if ("fine" %in% input$overlay) {
                    message("meshplot.expr: fine overlay")
                    plot(fine(), add=TRUE)
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

        limits <- reactiveValues(input.xlim=c(0,1),
                                 input.ylim=c(0,1),
                                 inputplot.xlim=c(0,1),
                                 inputplot.ylim=c(0,1))
        observe({
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
            limits$input.xlim <- r

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
            limits$input.ylim <- r
        })
        observe({
            limits$inputplot.xlim <- limits$input.xlim + c(-1,1)*input$offset[2]
            limits$inputplot.ylim <- limits$input.ylim + c(-1,1)*input$offset[2]
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
            plot(NA, type="n", main="",
                 xlim=limits$inputplot.xlim, ylim=limits$inputplot.ylim,
                 xlab="Easting", ylab="Northing",
                 asp=1)
            points(mesh.loc(), col=2)
            points(boundary.loc(), col=4)
            bnd <- boundary()
            for (k in 1:2) {
                if (inherits(bnd[[k]], "inla.mesh.segment")) {
                    lines(boundary()[[k]], col=k, add=TRUE)
                } else if (inherits(bnd[[k]], "list")) {
                    lapply(boundary()[[k]], function(x) lines(x, col=k, add=TRUE))
                }
            }
        })

        output$inputplot <-
            renderPlot(inputplot.expr, quoted=TRUE, height=800, units="px")

        output$meshmeta <- renderTable({meshmeta()}, digits=3, rownames=TRUE)

        clicks <- reactiveValues(meshplot.loc=NULL)
        observeEvent(input$meshplot.click, {
            message("meshplot.loc")
            if (is.null(input$meshplot.click)) {
                clicks$meshplot.loc <- NULL
            } else {
                clicks$meshplot.loc <- cbind(input$meshplot.click$x, input$meshplot.click$y)
            }
            message(paste("meshplot.loc: meshplot.loc =", paste(clicks$meshplot.loc, collapse=", ")))
        })
        meshplot.A.mesh <- reactive({
            message("meshplot.A.mesh")
            A <- NULL
            if (input$assess && !is.null(clicks$meshplot.loc) && !is.null(mesh())) {
                message(paste("meshplot.A.mesh:",
                              paste0("meshplot.loc = ", paste(clicks$meshplot.loc, collapse=", "))))
                A <- INLA::inla.spde.make.A(mesh(), clicks$meshplot.loc)
                if (!(sum(A) > 0)) {
                    A <- NULL
                }
            }
            A
        })
        meshplot.A.fine <- reactive({
            message("meshplot.A.fine")
            A <- NULL
            if (input$assess && !is.null(clicks$meshplot.loc) && !is.null(fine())) {
                message(paste("meshplot.A.fine:",
                              paste0("meshplot.loc = ", paste(clicks$meshplot.loc, collapse=", "))))
                A <- INLA::inla.spde.make.A(fine(), clicks$meshplot.loc)
                if (!(sum(A) > 0)) {
                    A <- NULL
                }
            }
            A
        })

        click.information <- reactive({
            message("click.information")
            if (input$assess.advanced) {
                if (is.null(clicks$meshplot.loc) ||
                    is.null(meshplot.A.mesh()) ||
                    is.null(meshplot.A.fine())) {
                    out <- rbind(cbind(NA, NA, NA), cbind(NA, NA, NA), cbind(NA, NA, NA))
                } else {
                    message("click.information!")
                    A.mesh <- meshplot.A.mesh()
                    A.fine <- meshplot.A.fine()
                    sd.bound.mesh <- as.vector(A.mesh %*% diag(mesh.S())^0.5)
                    sd.bound.fine <- as.vector(A.fine %*% diag(fine.S())^0.5)
                    sd.mesh <- Matrix::rowSums(A.mesh * (A.mesh %*% mesh.S()))^0.5
                    sd.fine <- Matrix::rowSums(A.fine * (A.fine %*% fine.S()))^0.5
                    out <- rbind(cbind(sd.mesh, sd.fine, sd.mesh/sd.fine),
                                 cbind(sd.bound.mesh,
                                       sd.bound.fine,
                                       sd.bound.mesh/sd.bound.fine),
                                 cbind(sd.mesh/sd.bound.mesh,
                                       sd.fine/sd.bound.fine,
                                       sd.mesh/sd.bound.mesh / (sd.fine/sd.bound.fine)))
                }
                rownames(out) <- c("SD", "SD bound", "SD/bound")
                colnames(out) <- c("Mesh", "Fine", "Mesh/Fine")
                out <- as.data.frame(out)
            } else {
                if (is.null(clicks$meshplot.loc) ||
                    is.null(meshplot.A.mesh())) {
                    out <- rbind(NA, NA, NA)
                } else {
                    message("click.information!")
                    A.mesh <- meshplot.A.mesh()
                    sd.bound.mesh <- as.vector(A.mesh %*% diag(mesh.S())^0.5)
                    sd.mesh <- Matrix::rowSums(A.mesh * (A.mesh %*% mesh.S()))^0.5
                    out <- rbind(sd.mesh, sd.bound.mesh, sd.mesh/sd.bound.mesh)
                }
                rownames(out) <- c("SD", "SD bound", "SD/bound")
                colnames(out) <- c("Mesh")
                out <- as.data.frame(out)
            }
            out
        })
        SD.information <- reactive({
            message("SD.information")
            if (!input$assess) {
                NULL
            } else {
                message("SD.information!")
                if (input$assess.advanced) {
                    out <- matrix(NA, 3,5)
                } else {
                    out <- matrix(NA, 1,5)
                }
                val <- as.vector(mesh.sd())
                if (!is.null(val)) {
                    out[1,1] <- min(val, na.rm=TRUE)
                    out[1,2:4] <- quantile(val, c(0.25, 0.5, 0.75), na.rm=TRUE)
                    out[1,5] <- max(val, na.rm=TRUE)
                    if (input$assess.advanced) {
                        val <- as.vector(fine.sd())
                        out[2,1] <- min(val, na.rm=TRUE)
                        out[2,2:4] <- quantile(val, c(0.25, 0.5, 0.75), na.rm=TRUE)
                        out[2,5] <- max(val, na.rm=TRUE)
                        val <- as.vector(mesh.sd() / fine.sd())
                        out[3,1] <- min(val, na.rm=TRUE)
                        out[3,2:4] <- quantile(val, c(0.25, 0.5, 0.75), na.rm=TRUE)
                        out[3,5] <- max(val, na.rm=TRUE)
                    }
                }
                rnames <- c("SD on Mesh",
                            "SD on Fine",
                            "SD relative deviation")
                if (input$assess.advanced) {
                    rownames(out) <- rnames
                } else {
                    rownames(out) <- rnames[1]
                }
                colnames(out) <- c("min", "25%", "50%", "75%", "max")
                out <- as.data.frame(out)
                out
            }
        })
        corr.information <- reactive({
            message("corr.information")
            if (!input$assess) {
                NULL
            } else {
                message("corr.information!")
                if (input$assess.advanced) {
                    out <- matrix(NA, 3,5)
                } else {
                    out <- matrix(NA, 1,5)
                }
                val <- as.vector(mesh.corr())
                if (!is.null(val)) {
                    out[1,1] <- min(val, na.rm=TRUE)
                    out[1,2:4] <- quantile(val, c(0.25, 0.5, 0.75), na.rm=TRUE)
                    out[1,5] <- max(val, na.rm=TRUE)
                    if (input$assess.advanced) {
                        val <- as.vector(fine.corr())
                        out[2,1] <- min(val, na.rm=TRUE)
                        out[2,2:4] <- quantile(val, c(0.25, 0.5, 0.75), na.rm=TRUE)
                        out[2,5] <- max(val, na.rm=TRUE)
                        val <- as.vector(mesh.corr() - fine.corr())
                        out[3,1] <- min(val, na.rm=TRUE)
                        out[3,2:4] <- quantile(val, c(0.25, 0.5, 0.75), na.rm=TRUE)
                        out[3,5] <- max(val, na.rm=TRUE)
                    }
                }
                rnames <- c("Correlation on Mesh",
                            "Correlation on Fine",
                            "Correlation deviation")
                if (input$assess.advanced) {
                    rownames(out) <- rnames
                } else {
                    rownames(out) <- rnames[1]
                }
                colnames(out) <- c("min", "25%", "50%", "75%", "max")
                out <- as.data.frame(out)
                out
            }
        })
        
        output$click.information <- renderTable({
            if (input$assess) {
                click.information()
            }
        }, digits=3, rownames=TRUE)
        output$field.information <- renderTable({
            if (input$assess.quantity == "corr") {
                corr.information()
            } else {
                SD.information()
            }
        }, digits=3, rownames=TRUE)

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
                              attr(fine(), "code"), "\n")
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

        output$input.loc.coordinates <- renderText({
            xy_str <- function(e) {
                if (is.null(e)) {
                    "NULL\n"
                } else {
                    paste0("x=", e$x, ", y=", e$y, "\n")
                }
            }
            paste0(xy_str(input$inputplot.hover))
        })
        output$meshplot.loc.coordinates <- renderText({
            xy_str <- function(e) {
                if (is.null(e)) {
                    "NULL\n"
                } else {
                    paste0("x=", e[1,1], ", y=", e[1,2], "\n")
                }
            }
            paste0(xy_str(clicks$meshplot.loc))
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
