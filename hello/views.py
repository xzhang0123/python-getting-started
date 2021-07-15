from django.shortcuts import render
from django.http import HttpResponse
import requests
from .models import Greeting
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import least_squares

from plotly.offline import plot
from plotly.graph_objs import Scatter
# Create your views here.
def index(request):
    #r = requests.get('http://httpbin.org/status/418')
    #print(r.text)
    #return HttpResponse('<pre>' + r.text + '</pre>')
    # return HttpResponse('Hello from Python!')
    return render(request, "index.html")


def db(request):

    greeting = Greeting()
    greeting.save()

    greetings = Greeting.objects.all()

    return render(request, "db.html", {"greetings": greetings})

def signupuser(request):
    # return render(request,"index.html")
    if request.method == "GET":
        return render(request,"signupuser.html",{'form':UserCreationForm()})
    elif request.method =="POST":
        if request.POST['password1'] ==request.POST['password2']:
            user = User.objects.create_user(request.POST['username'],password=request.POST['password1'])
            user.save()
        else:
            print("Password not consistent")

# all_data =pd.read_excel(r'C:\Users\zxiao\Google Drive\Xiao_Model\Experiments\R123_2020_Kidney\PTI_1ADP_all_data.xlsx')
# Int =all_data['C_GM_Avg'].dropna()
Int = [1, 0.9,0.7,0.65, 0.65,0.65,0.65,0.61,0.61,0.60,0.60,0.73,0.70,0.70,0.64,0.64,0.88]
F =0.096484
R =8.314e-3
T = 310.15
C0 = 1e-6
Ve = 1000;
Vx = 1;

beta = 0.33
alfa = 4.49
k0 = 8.15e4
def Con2Int(C):
    k=1
    k0 = 8.15e4
    k1 = 5.16e4
    k2 = 1.01e6
    k3 = 3.91e12
    I = k0*C/(1+k1*C+k2*C**2+k3*C**3)
    return I

def equations(vars,I):
    Ce,dPsi = vars
    #----------------
    Cx = Ce*np.exp(F*dPsi/(R*T))
    I0 = Con2Int(C0)
    Ie = Con2Int(Ce)
    Ix = Con2Int(Cx)
    eq1 =(Ve*Ie + Vx*Ix + beta*k0*alfa *(Ce+Cx)/2 )/(Ve*I0) -I
    eq2 = Ve*Ce + Vx*Cx +alfa *(Ce+Cx)/2 - Ve*C0
    return [eq1,10000*eq2]

def cal(request):
    y_data=[]
    for i in Int:
        I =i
        res = least_squares(equations, (1e-8,160), bounds = ((0, 0), (1, 240)),args=(I,))
    #     print(res.x)
        y_data.append(res.x[1])

    # y_data = [x**2 for x in x_data]
    plot_div = plot([Scatter(x=list(range(len(Int))), y=y_data,
                        # title='Converted Mito-Membrane Potential',
                        mode='lines', name='test',
                        opacity=0.8, marker_color='green')],
               output_type='div')
    plot_int_div = plot([Scatter(x=list(range(len(Int))), y=Int,
                        # title='Dye intensity data',
                        mode='lines', name='test',
                        opacity=0.8, marker_color='green')],
               output_type='div')
    return render(request, "cal.html", context={'plot_int_div': plot_int_div,'plot_div': plot_div})
