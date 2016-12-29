import os

folder = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(folder, '../VERSION')) as f:
    version = f.read().strip()

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
]
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = 'gig'
copyright = '2016, Danilo Horta'
author = 'Danilo Horta'
release = version
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
html_theme = 'default'
html_static_path = ['_static']
htmlhelp_basename = 'gigdoc'
latex_elements = {}
latex_documents = [
    (master_doc, 'gig.tex', 'gig Documentation',
     'Danilo Horta', 'manual'),
]
man_pages = [
    (master_doc, 'gig', 'gig Documentation',
     [author], 1)
]
texinfo_documents = [
    (master_doc, 'gig', 'gig Documentation',
     author, 'gig', 'One line description of project.',
     'Miscellaneous'),
]
