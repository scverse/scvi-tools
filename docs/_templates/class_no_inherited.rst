{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. add toctree option to make autodoc generate the pages

.. autoclass:: {{ objname }}

   {% block methods %}
   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
      :toctree: .
   {% for item in methods %}
      {%- if item != '__init__' and item not in inherited_members%}
        ~{{ fullname }}.{{ item }}
      {%- endif -%}

   {%- endfor %}
   {% endif %}
   {% endblock %}
