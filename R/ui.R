


ui <- fluidPage(

  shinyjs::useShinyjs(),

  # include custom design
  includeCSS("R/www/style.css"),

  # Application title
  #titlePanel("HD prior"),

  h2(""),

  fluidRow(

    # option panels
    column(4,
           shinyBS::bsCollapse(id = "sidebarpanel", multiple = TRUE, open = c("choose_structure", "choose_priors"),#, "closing"),
                      shinyBS::bsCollapsePanel(
                        #style = "danger",
                        title = "Modify tree structure",
                        value = "choose_structure",

                        column(12,
                        disabled(
                          detach_button(),
                          attach_button()
                        ),
                        disabled(
                          merge_button(),
                          remove_button()
                        ))
                      ), # end collapsable panel

                      shinyBS::bsCollapsePanel(
                        title = "Change priors",
                        value = "choose_priors",

                        div(#style = "min-height: 122px; max-height = 122px",
                        fluidRow(
                          column(12,
                                 h5("Prior on split"),
                        column(6, style = 'padding: 0px;',
                        disabled(
                          radioButtons("weight_prior",
                                       NULL,
                                       choices = list("Dirichlet" = "dirichlet", "PC prior" = "pc"),
                                       selected = "dirichlet")
                        )
                        ),
                        column(6, style = 'padding: 0px;',
                               div(class = "basemodelradio",
                                   radioButtons("basemodel", "Base model",
                                                choiceValues = c("mod1", "mod2", "comb"),
                                                choiceNames = c("100% A + 0% B", "0% A + 100% B", "25% A + 75% B"),
                                                selected = "mod1"
                                   )),
                               numericInput("median",
                                            "Median",
                                            value = 0.5, min = 0, max = 1, step = 0.05),
                               numericInput("conc_param",
                                            "Concentration",
                                            value = 0.5, min = 0.5, max = 1, step = 0.05)
                        )))
                        ),

                        h4(""),
                        fluidRow(
                          column(12,
                                 h5(textOutput("variance_prior_name")),
                        column(6, style = 'padding: 0px;',
                        disabled(
                          radioButtons("var_prior",
                                       NULL,
                                       choices = list("Jeffreys'" = "jeffreys",
                                                      "PC prior" = "pc0",
                                                      "Inverse gamma" = "invgam",
                                                      "Half-Cauchy (st.dev.)" = "hc",
                                                      "Half-normal (st.dev.)" = "hn"),
                                       selected = "jeffreys")
                        )
                        ),
                        column(6, style = 'padding: 0px;',
                          hidden(numericInput("par1", "", value = 1, min = 0, max = Inf, step = 1)),
                          hidden(numericInput("par2", "", value = 1, min = 0, max = 1, step = 0.05))
                        )
                          )
                      )), # end collapsable panel

                      shinyBS::bsCollapsePanel(
                        title = "Parameterization",
                        value = "param_expressions",
                        div(style = "overflow-x: scroll", uiOutput("param_eq")),
                        checkboxInput("param_type_eq", label = "Function of variances", value = TRUE)
                      ), # end collapsable panel

                      shinyBS::bsCollapsePanel(
                        title = "Priors",
                        value = "prior_distributions",
                        div(style = "overflow-x: scroll", uiOutput("prior_distributions"))
                      ) # end collapsable panel
           ),
           fluidPage(
             theme = "default",
             uiOutput("begin_quit_guide"),
             actionButton("reset", "Reset", style = "width:70px", class = "button1", title = "Reset everything to the settings you started with."),
             actionButton("close", "Close", style = "width:70px", class = "button1", title = "Close application. The application stores your settings when closing. You can also just close the window, this will still save your settings.")
             # actionButton("close_and_run", "Close and run", style = "width:140px", class = "button1", title = "Close application and start inference with INLA immediately.")
           )

    ), # end sidebarpanel


    # main panel
    mainPanel(width = 8,

              h4(
                uiOutput("model_eq"),
              ),

              fluidRow(
                column(12,
                       visNetwork::visNetworkOutput("graph"),
                       div(class = "fix_height_12", textOutput("chosen_node")),
                       h5(""),
                       textOutput("intrinsic"),
                       h5("")
                ),

                # plots
                column(width = 12, # style = "border: 1px solid black; background-color: white",
                       uiOutput("plot_or_next_step_button"),

                       textOutput("no_pc"),
                       
                       textOutput("guide_message"),

                       disabled(plotOutput("all_plots"))

                ) # end column

              )
    ) # end mainpanel
  ) # end fluidrow

)




