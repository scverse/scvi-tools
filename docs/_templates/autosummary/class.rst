{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. add toctree option to make autodoc generate the pages

.. autoclass:: {{ objname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
   {% for item in attributes %}
      ~{{ fullname }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}
   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
   {% for item in methods %}
      {%- if item != '__init__' %}
        ~{{ fullname }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes_documentation %}
   {% if attributes %}

   {% for item in attributes %}
   .. autoattribute:: {{ item }}
   {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block methods_documentation %}
   {% if methods %}

   {% for item in methods %}
   .. automethod:: {{ item }}
      {%- if item != '__init__' %}
        .. automethod:: {{ item }}
      {%- endif -%}
   {%- endfor %}

   {% endif %}
   {% endblock %}
