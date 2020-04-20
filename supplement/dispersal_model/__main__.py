"""Run the dispersal model."""

import argparse

import dispersal_model.dispersal as d


# # Appendix: Running the Model
# The model has a command line interface to set most parameters and run it.
# Running sets the global `params` variable. Then the model is initialized
# according to the parameters, and ran according to the schedule for
# 15000 simulated years or until all families have died.
def main() -> None:
    """Command line interface entry point."""
    parser = argparse.ArgumentParser(
        description=__doc__)
    parameter_settings = parser.add_argument_group("Model Parameters")
    for parameter, field in d.ParameterSetting.__dataclass_fields__.items():
        parameter_settings.add_argument(
            "--{:s}".format(parameter.replace("_", "-")),
            type=field.type,
            default=field.default
        )

    parser.add_argument(
        "--log", type=argparse.FileType('w'), default='log',
        help="Logfile for observations")
    parser.add_argument(
        "--report_resources", action="store_true", default=False,
        help="Include patch resources in log file")

    args = parser.parse_args()
    d.params = d.ParameterSetting(**{
        name: value
        for name, value in vars(args).items()
        if name in d.ParameterSetting.__dataclass_fields__})

    s = d.initialization()

    for i in range(15000 * 2):
        s = d.step(s)
        d.observation(s, args.log, report_resources=args.report_resources)


if __name__ == "__main__":
    main()
