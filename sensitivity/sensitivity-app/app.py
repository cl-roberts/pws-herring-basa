from shiny import App, ui, render, reactive
import os
import atexit

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
    ui.layout_columns(
        ui.value_box(
            "Interesting bit", ui.output_ui("interesting_bit"),
        ),
        ui.value_box(
            "Fascinating result", ui.output_ui("fascinating_result"), 
        ),
        ui.value_box(
            "Fabulous insight",
            ui.output_ui("fabulous_insight"),
        ),
        fill=False,
    ),
    ui.layout_columns(
        ui.card(
            ui.card_header("Comparison to base model biomass"),
            ui.output_image("compare_biomass"),
            full_screen=True,
        ),
    ),
)

def server(input, output, session):

    @reactive.effect
    @reactive.event(input.action_button)
    def run_sensitivity():
        sensitivity.modify_mortality(new_value=input.m(), dir_model=dir_model)
        return_value = sensitivity.run_basa(dir_sensitivity=dir_sensitivity)
        sensitivity.plot_sensitivity(dir_outputs, dir_outputs + "/sensitivity_plot.png")
        return return_value

    @render.text
    def interesting_bit():
        print(dir_model)
        return input.m()

    @reactive.file_reader(dir_outputs + "/outputs-for-management.csv")
    @render.image
    def compare_biomass():
        dir = dir_outputs + "/sensitivity_plot.png"

        if os.path.exists(dir):
            return {"src": dir, "width": "500px"}
        else:
            return 

app = App(app_ui, server)