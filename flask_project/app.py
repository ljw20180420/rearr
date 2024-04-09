from flask import Flask

# change Jinja2 delimiter to solve the comflicts with vue3
class VueFlask(Flask):
    jinja_options = Flask.jinja_options.copy()
    jinja_options.update(dict(
        block_start_string='{b{s',
        block_end_string='b}e}',
        variable_start_string='{v{s',
        variable_end_string='v}e}',
        comment_start_string='{c{s',
        comment_end_string='c}e}',
    ))
