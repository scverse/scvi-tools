{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. add toctree option to make autodoc generate the pages

.. autoclass:: {{ objname }}
   :show-inheritance:

{% block attributes %}
{% if attributes %}
Attributes
~~~~~~~~~~

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
~~~~~~~~

.. autosummary::
{% for item in methods %}
    {%- if item != '__init__' and item not in inherited_members%}
    ~{{ fullname }}.{{ item }}
    {%- endif -%}

{%- endfor %}
{% endif %}
{% endblock %}

{% block attributes_documentation %}
{% if attributes %}

{% for item in attributes %}
{{ item }}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoattribute:: {{ [objname, item] | join(".") }}
{%- endfor %}

{% endif %}
{% endblock %}

{% block methods_documentation %}
{% if methods %}

{% for item in methods %}
{%- if item != '__init__' and item not in inherited_members%}
{{ item }}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: {{ [objname, item] | join(".") }}
{%- endif -%}
{%- endfor %}

{% endif %}
{% endblock %}
