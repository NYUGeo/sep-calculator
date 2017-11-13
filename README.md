# Seismic Earth Pressure Calculator

### Click [here](http://cue3.engineering.nyu.edu:5010) for the live web application, or read below to run locally.

<br>

[![One layer](media/l1-screenshot.png)](http://cue3.engineering.nyu.edu:5010)

If you find our work useful, please cite us:

> Machairas N, Iskander M, Omidvar M (2018) “Interactive Web Application for Computing Seismic Earth Pressure” *Proceedings GEESD V 2018*, ASCE (in review)


<br>

## What is the Seismic Earth Pressure (SEP) Calculator

An interactive web application has been developed to compute seismic earth pressure. The application employs an expanded form of Rankine's classic earth pressure solution. The application computes seismic active earth pressure behind rigid walls supporting c–φ backfill considering both wall inclination and backfill slope. The lateral earth pressure formulation is based on Rankine's conjugate stress concept. The developed expression can be used for the static and pseudo-static seismic analyses of c–φ backfill. The results based on the proposed formulations are found to be identical to those computed with the Mononobe–Okabe method for cohesionless soils, provided the same wall friction angle is employed. For c–φ soils, the formulation yields comparable results to available solutions for cases where a comparison is feasible. The application eliminates the need for design charts, since it can accommodate any set of user-defined kinematically admissible design parameters.

<br>

## Run on your computer

#### Requirements

All programming was done in Python 3.6. You may run the program in any operating system with Python 3.4 or above installed. Required packages for full functionality: `numpy` and `bokeh`. Use `git clone` (recommended) or download ZIP to get this repository on your computer. From within the local repo directory, run the following in your command prompt (terminal) to install the necessary packages with compatible versions.

```
pip install -r requirements.txt
```

There are two ways to perform calculations:
- A) manually in the command prompt (terminal) and/or
- B) interactively in the web application.


### A) Run in Command Prompt (w/ example)

Development of the core algorithms of the SEP calculator is closely following:

>Iskander, M., Chen, Z., Omidvar, M., Guzman, I., and Elsherif, O. (2013). “Active static and seismic earth pressure for c–φ soils.” *Soils and Foundations*, 53(5), 639–652.

The following example will replicate Figure 5 of the aforementioned journal article.

Launch a Python shell from within the project folder. Import the core SEP class with:

```python
>>> from sep_core import sep
```

**Input parameters:** The current version is working with S.I. units only. The input parameters for the example shown in Figure 5 are:

- kh = 0.2
- kv = -0.1
- &#969; = 20&#176;
- &#946; = 15&#176;
- &#966; = 30&#176;
- &#947; = 23 kN/m&#179;
- c = 20 kPa
- H = 15 m

Instantiate the Python SEP class with:

```python
>>> figure5 = sep(0.2, -0.1, 20, 15, 30, 23, 20, 15)
```

All necessary equations have been programmed as class methods and the values are returned in terms of depth from top of the wall, *Zw*. Therefore, at depth *Zw* = 9 meters, the following methods replicate the table in Figure 5:

```python
>>> figure5.zl(9)
9.577599952283208

>>> figure5.z(9)
9.8777297730653579

>>> figure5.Ja(9)
157.46389426135741

>>> figure5.alpha_a(9, degrees=True)
32.234474051287712

>>> figure5.Ka(9)
0.79140861884533187

>>> figure5.sigma_a(9)
179.79837097166856

>>> figure5.sigma_AEH(9)
110.11418922580418
```


### B) Run interactively

From within the `apps` folder, for one layer, launch the web application with:

```
bokeh serve --show sep-calculator-l1.py
```

And for two layers:

```
bokeh serve --show sep-calculator-l2.py
```
