<h1>Workflow for discovering, querying and visualizing astronomy data in the VO</h1>


```python
from astropy.coordinates import SkyCoord
from hips import WCSGeometry, make_sky_image
from hips import HipsSurveyProperties
from ipywidgets import Layout, Box, widgets

```

<hr>

# Data Discovery (VO)


```python
import pyvo
```


```python
from pyvo.registry import search as regsearch
```

## Find all TAP Services with 'quasars'


```python
services = regsearch(keywords=['quasar'], servicetype='tap')
```


```python
print (services)
```

## Find all TAP Services with keyword "ukidss"


```python
services = regsearch(keywords=['ukidss'], servicetype='tap')
print (services)
```


```python
ivoid=b'ivo://wfau.roe.ac.uk/ukidssdr9-dsa'
matching_row_indices = [i for i, val in enumerate(services["ivoid"]==ivoid) if val]
first_row_index = matching_row_indices[0] if len(matching_row_indices) > 0 else None
if (first_row_index==None):
    print ("No Matching services")
else:
    record = services.getrecord(first_row_index)
tap_url = record["access_url"].decode('UTF-8')
```


```python
print("TAP URL: %s \n" % tap_url)
```

<hr>

# Data Access (TAP, Astropy)


```python
from pyvo.dal import tap
service = tap.TAPService(tap_url)

query_text = """
SELECT TOP 500 
sourceID, ra,dec FROM 
    lasSource
ORDER by dec
"""


from astroquery.utils.tap.core import TapPlus
service = TapPlus(url=tap_url)
job = service.launch_job(query_text)
table = job.get_results()
table
```

<hr>

# Image Access

## HIPS Image


```python
url = 'http://surveys.roe.ac.uk/hips71/LAS/LAS_Y_OUT//properties'
hips_survey = 'CDS/P/DSS2/red'
```


```python
print (url)
```


```python
hips_survey_property = HipsSurveyProperties.fetch(url)
```


```python
geometry = WCSGeometry.create(skydir=SkyCoord(217.9, -2.3, unit='deg', frame='fk5'),width=180, height=180, fov='0.02 deg',coordsys='icrs', projection='TAN')
from hips import make_sky_image
result = make_sky_image(geometry, hips_survey_property, 'fits') 
```


```python
result.plot()
```


```python
#result.image
```


```python
import matplotlib.pyplot as plt
```


```python
from astropy.visualization.mpl_normalize import simple_norm
```


```python
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
norm = simple_norm(result.image, 'sqrt', min_percent=1, max_percent=99)
im = ax.imshow(result.image, origin='lower', norm=norm, cmap='gray_r')
fig.colorbar(im)
```

<hr>

# Visualization

## Scatter Plot using Matplotlib


```python
%matplotlib notebook
from ipywidgets import *
import numpy as np
import matplotlib.pyplot as plt

x = table["ra"]
y = table["dec"]
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
line = ax.scatter(x, y)

plt.xlabel("$ra$", fontsize=20)
plt.ylabel("$decl$", fontsize=20)

def update(w = 1.0):
    fig.canvas.draw()

interact(update);
```


```python
from plotly import __version__
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
from plotly.graph_objs import *
import plotly.graph_objs as go

import numpy as np

init_notebook_mode() # Important
```

## Scatter Plot using Plotly


```python

N = 100000
trace = go.Scattergl(
    x = table["ra"],
    y = table["dec"],
    mode = 'markers',
    marker = dict(
        color = '#FFBAD2',
        line = dict(width = 1)
    )
)
data = [trace]
iplot(data, filename='compare_webgl')
```

## 3d Plot using Plotly


```python
x, y, z = np.random.multivariate_normal(np.array([0,0,0]), np.eye(3), 200).transpose()
trace1 = go.Scatter3d(
    x=table["ra"],
    y=table["dec"],
    z=z,
    mode='markers',
    marker=dict(
        size=12,
        line=dict(
            color='rgba(217, 217, 217, 0.14)',
            width=0.5
        ),
        opacity=0.8
    )
)

data = [trace1]
layout = go.Layout(
    margin=dict(
        l=0,
        r=0,
        b=0,
        t=0
    )
)
fig = go.Figure(data=data, layout=layout)
iplot(fig, filename='simple-3d-scatter')
```

<hr>

## Aladin Lite


```python
import ipyaladin as ipyal
```


```python
aladin=ipyal.Aladin(target='12 07 57.417 -03 58 16.02', fov=14, survey='http://surveys.roe.ac.uk/hips71/LAS/LAS_Y_OUT//properties')
aladin
```


```python
aladin.add_table(table)

```

## Linked Aladin Lite Views



```python
from ipywidgets import Layout, Box, widgets
a = ipyal.Aladin(layout=Layout(width='33.33%'), target='09 14 42.373 +00 29 43.26', fov=0.3)
b = ipyal.Aladin(layout=Layout(width='33.33%'), survey='P/DSS2/red')
c = ipyal.Aladin(layout=Layout(width='33.33%'), survey='P/2MASS/color')

# synchronize target between 3 widgets
widgets.jslink((a, 'target'), (b, 'target'))
widgets.jslink((b, 'target'), (c, 'target'))

# synchronize FoV (zoom level) between 3 widgets
widgets.jslink((a, 'fov'), (b, 'fov'))
widgets.jslink((b, 'fov'), (c, 'fov'))

items = [a, b, c]

box_layout = Layout(display='flex',
                    flex_flow='row',
                    align_items='stretch',
                    border='solid',
                    width='100%')
box = Box(children=items, layout=box_layout)
box
```

## Aladin Lite On-click action



```python
aladin = ipyal.Aladin(layout=Layout(width='70%'), target='M 36', fov=0.3)
info = widgets.HTML()


box_layout = Layout(display='flex',
                    flex_flow='row',
                    align_items='stretch',
                    width='100%')
box = Box(children=[aladin, info], layout=box_layout)
box
```


```python

import requests
def process_result(data):
    info.value = ''
    ra = data['ra']
    dec = data['dec']
    radius = min(aladin.fov / 10, 5)
    query = """SELECT TOP 1 main_id, ra, dec, DISTANCE(POINT('ICRS', %f, %f), POINT('ICRS', ra, dec)) as d FROM basic
               WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', %f, %f, %f))=1
               ORDER BY d ASC""" % (ra, dec, ra, dec, aladin.fov / 10)
    
    r = requests.get('http://simbad.u-strasbg.fr/simbad/sim-tap/sync', params={'query': query, 'request': 'doQuery', 'lang': 'adql', 'format': 'json', 'phase': 'run'})
    obj_name = ''
    obj_coo = None
    obj_data = r.json()['data']
    if len(obj_data)==0:
        return
    
    obj_name = obj_data[0][0]
    obj_coo = [obj_data[0][1], obj_data[0][2]]
    sed_img = '<img src="http://cdsportal.u-strasbg.fr/cgi-bin/PhotVizPreview/plot?ra=%f&dec=%f&radius_arcsec=5&w=200&h=150&point_size=4">' % (obj_coo[0], obj_coo[1])
    info.value =  '<h2>%s</h2><br><br>%s' % (obj_name, sed_img)
    
aladin.add_listener('click', process_result)

```

## Aladin Lite Selection



```python

from ipyaladin import Aladin
from ipywidgets import Layout, Box, widgets

aladin = Aladin(layout=Layout(width='50%'), target='12 06 41.797 -03 47 54.95', fov=3.26,survey='http://surveys.roe.ac.uk/hips71/LAS/LAS_Y_OUT//properties')


button = widgets.Button(description="Select")
def on_button_clicked(b):
    aladin.rectangular_selection()

button.on_click(on_button_clicked)
table_info = widgets.HTML(layout=Layout(height='auto', overflow='auto'))


box_layout = Layout(display='flex',
                    flex_flow='row',
                    align_items='stretch',
                    width='100%',
                    position='relative',
                    overflow='hidden',
                    height='100vh',
                    margin='-100px 0 0 0',
                    padding='100px 0 0 0 '
                   )
box = Box(children=[aladin, button, table_info], layout=box_layout)
box
```


```python
from astroquery.simbad import Simbad
import astropy.units as u

aladin.add_table(table)

def process_result(sources):
    s = '<table border="1">'
    s += '<tr><th>Source ID</th><th>RA</th><th>DEC</th></tr>'
    for source in sources:
        s += '<tr><td>%s</td><td>%s</td><td>%s</td></tr>' % (source['data']['sourceID'], source['data']['ra'], source['data']['dec'])
    s += '</table>'
    table_info.value = s
    
aladin.add_listener('select', process_result)
```


```python

```
