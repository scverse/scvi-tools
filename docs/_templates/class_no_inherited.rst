{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. add toctree option to make autodoc generate the pages

.. autoclass:: {{ objname }}
   :show-inheritance:

   {% block attributes %}
   {% if attributes %}
   Attributes
   ^^^^^^^^^^

   .. autosummary::
    {% for item in attributes %}
        {%- if item not in inherited_members%}
            ~{{ fullname }}.{{ item }}
        {%- endif -%}
    {%- endfor %}
    {% endif %}
    {% endblock %}


   {% block methods %}
   {% if methods %}
   Methods
   ^^^^^^^

   .. autosummary::
   {% for item in methods %}
      {%- if item != '__init__' and item not in inherited_members%}
        ~{{ fullname }}.{{ item }}
      {%- endif -%}

   {%- endfor %}
   {% endif %}
   {% endblock %}
