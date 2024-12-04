from shiny import App, ui, render, reactive
import os

from sensitivity.src import sensitivity

dir_base = os.getcwd()
dir_sensitivity = dir_base + "/sensitivity"
dir_model = dir_sensitivity + "/model"
dir_outputs = dir_sensitivity + '/data_outputs'

app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_slider(
            "m",
            "Natural Mortality",
            min=0.0,
            max=0.5,
            value=0.25,
        ),
        ui.input_action_button("action_button", "Run Model"),
    ),
    ui.card(
        ui.card_header("Welcome"),
        ui.output_ui("message_box"),
    ),
    ui.card(
        ui.card_header("Comparison to base model biomass"),
        ui.output_image("compare_biomass"),
    ),
    ui.card(
        ui.card_header("Error analysis"),
        ui.output_image("error_analysis"),
        ui.layout_columns(
            ui.value_box(
                "Mean Percent Error", 
                ui.output_ui("mean_error"),
            ),
            ui.value_box(
                "Min Percent Error", 
                ui.output_ui("min_error"), 
            ),
            ui.value_box(
                "Max Percent Error",
                ui.output_ui("max_error"),
            ),
        ),
    ),
)


def server(input, output, session):

    # return_value = reactive.value(0)

    @render.text 
    def message_box():
        if os.path.exists(dir_outputs + "/error_plot.png") == False:
            return "Select a value for natural mortality and run model to get started"
        # if return_value() == 1:
        #     return "Model failed to converge. Try a different value for natural mortality"


    @reactive.effect
    @reactive.event(input.action_button)
    def run_sensitivity():
        sensitivity.modify_mortality(new_value=input.m(), dir_model=dir_model)
        sensitivity.run_basa(dir_sensitivity=dir_sensitivity)
        # value = sensitivity.run_basa(dir_sensitivity=dir_sensitivity)
        # return_value.set(value)
        sensitivity.plot_sensitivity(dir_outputs, dir_outputs + "/sensitivity_plot.png")
        sensitivity.sens_comparison(dir_outputs)
        sensitivity.plot_sens_comparison(dir_outputs, dir_outputs + "/error_plot.png")
        return 

    @reactive.file_reader(dir_outputs + "/outputs-for-management.csv")
    @render.text
    def mean_error():
        if os.path.exists(dir_outputs + "/error_plot.png"):
            mean_error = sensitivity.sens_comparison(dir_outputs).at[0, 'Mean Perc. Error']
            return str(round(mean_error, 2)) + "%"
        else:
            return

    @reactive.file_reader(dir_outputs + "/outputs-for-management.csv")
    @render.text
    def min_error():
        if os.path.exists(dir_outputs + "/error_plot.png"):
            min_error = sensitivity.sens_comparison(dir_outputs).at[0, 'Min Perc. Error']
            return str(round(min_error, 2)) + "%"
        else:
            return
    
    @reactive.file_reader(dir_outputs + "/outputs-for-management.csv")
    @render.text
    def max_error():
        if os.path.exists(dir_outputs + "/error_plot.png"):
            max_error = sensitivity.sens_comparison(dir_outputs).at[0, 'Max Perc. Error']
            return str(round(max_error, 2)) + "%"
        else:
            return
    
    @reactive.file_reader(dir_outputs + "/outputs-for-management.csv")
    @render.image
    def compare_biomass():
        dir = dir_outputs + "/sensitivity_plot.png"

        if os.path.exists(dir):
            return {"src": dir, "width": "500px"}
        else:
            return 

    @reactive.file_reader(dir_outputs + "/outputs-for-management.csv")
    @render.image
    def error_analysis():
        dir = dir_outputs + "/error_plot.png"

        if os.path.exists(dir):
            return {"src": dir, "width": "500px"}
        else:
            return 

app = App(app_ui, server)