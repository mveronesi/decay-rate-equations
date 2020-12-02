import ipywidgets as widgets
from ipywidgets import interact, interactive, interactive_output
from ipywidgets import Layout, Button, Box, VBox, HBox
from plot_amplitudes import plot_amplitudes
from declare_acc import a_acc, n_acc, b_acc, beta_acc, cutoff_acc

# bmix - cp parameters
dm=widgets.BoundedFloatText(value=17.749, min=0, max=20, step=0.001, continuous_update=False, description=r'$\Delta m_{s}$ [ps$^{-1}$]')
dg=widgets.BoundedFloatText(value=0.083, min=0, max=1, step=0.001, continuous_update=False, description=r'$\Delta \Gamma_{s}$ [ps$^{-1}$]')
gs=widgets.BoundedFloatText(value=0.6608, min=0, max=1, step=0.001, continuous_update=False, description=r'$\Gamma_{s}$ [ps$^{-1}$]')
r=widgets.BoundedFloatText(value=0, min=0, max=1, step=0.1, continuous_update=False, description=r'$r$')
delta=widgets.BoundedFloatText(value=0, min=0, max=360, step=1, continuous_update=False, description=r'$\delta$ [$\circ$]')
gamma=widgets.BoundedFloatText(value=0, min=0, max=360, step=1, continuous_update=False, description=r'$\gamma$ [$\circ$]')
beta=widgets.BoundedFloatText(value=0, min=-1000, max=1000, step=1, continuous_update=False, description=r'$\beta_{s}$ [mrad]')
# afs, cpt
afs=widgets.BoundedFloatText(value=0, min=-1, max=1, step=0.001, continuous_update=False, description=r'$a_{fs}$')
rez=widgets.BoundedFloatText(value=0, min=-1, max=1, step=0.001, continuous_update=False, description=r'$Re(z)$')
imz=widgets.BoundedFloatText(value=0, min=-1, max=1, step=0.001, continuous_update=False, description=r'$Im(z)$')
# tagging
omega=widgets.BoundedFloatText(value=0, min=0, max=0.5, step=0.01, continuous_update=False, description=r'$\omega_{tag}$')
d_omega=widgets.BoundedFloatText(value=0, min=-0.5, max=0.5, step=0.001, continuous_update=False, description=r'$\Delta\omega_{tag}$')
eff=widgets.BoundedFloatText(value=1, min=0, max=1, step=0.01, continuous_update=False, description=r'$\varepsilon_{tag}$')
d_eff=widgets.BoundedFloatText(value=0, min=-0.5, max=0.5, step=0.001, continuous_update=False, description=r'$\Delta\varepsilon_{tag}$')
sigma_t=widgets.BoundedFloatText(value=0, min=0, max=1, step=0.001, continuous_update=False, description=r'$\sigma_{t}$ [ps]')
# asymmetries
a_prod=widgets.BoundedFloatText(value=0, min=-1, max=1, step=0.001, continuous_update=False, description=r'$a_{prod}$')
a_det=widgets.BoundedFloatText(value=0, min=-1, max=1, step=0.001, continuous_update=False, description=r'$a_{det}$')
# select
b_f = widgets.Checkbox(value=True,description=r'$B$$\to$$f$',disabled=False)
bbar_f = widgets.Checkbox(value=False,description=r'$\overline{B}$$\to$$f$',disabled=False)
bbar_fbar = widgets.Checkbox(value=False,description=r'$\overline{B}$$\to$$\overline{f}$',disabled=False)
b_fbar = widgets.Checkbox(value=False,description=r'$B$$\to$$\overline{f}$',disabled=False)
u_f = widgets.Checkbox(value=False,description=r'$U$$\to$$f$',disabled=False)
u_fbar = widgets.Checkbox(value=False,description=r'$U$$\to$$\overline{f}$',disabled=False)
sum_f = widgets.Checkbox(value=False,description=r'Sum $f$',disabled=False)
sum_fb = widgets.Checkbox(value=False,description=r'Sum $\bar{f}$',disabled=False)
# decay rate coefficients
xmin=widgets.BoundedFloatText(value=0, min=0, max=20, step=1, continuous_update=False, description=r'$t_{min}$ [ps]')
xmax=widgets.BoundedFloatText(value=5, min=0, max=20, step=1, continuous_update=False, description=r'$t_{max}$ [ps]')
y_cosh=widgets.BoundedFloatText(value=0.25, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{\cosh}$ [a.u.]')
y_sinh=widgets.BoundedFloatText(value=0.02, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{\sinh}$ [a.u.]')
y_cos=widgets.BoundedFloatText(value=0.25, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{\cos}$ [a.u.]')
y_sin=widgets.BoundedFloatText(value=0.25, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{\sin}$ [a.u.]')
y_dg=widgets.BoundedFloatText(value=0.5, min=0, max=10, step=0.1, continuous_update=False, description=r'$y_{d\Gamma/dt}$ [a.u.]')
name=widgets.Text(value='plots/coeffs.eps',continuous_update=False,description='Save as:')
save=widgets.ToggleButton(value=False,description='Save')
tagging =widgets.Dropdown(
    options=['MESON', 'TAGGER'],
    value='MESON',
    description='Tagging',
    disabled=False,
)
k_acc = widgets.Checkbox(value=False,description=r'Acceptance',disabled=False)
k_tag = widgets.Checkbox(value=False,description=r'Tagging',disabled=False)
## Decay Rate Coefficients
wbox_amp = HBox([VBox([dm,dg,gs,r,delta,gamma,beta,afs,rez,imz]),
                 VBox([k_tag,a_prod,a_det,omega,d_omega,eff,d_eff,tagging]),
                 VBox([k_acc,a_acc,n_acc,b_acc,beta_acc,cutoff_acc,sigma_t]),
                 VBox([xmin,xmax,y_cosh,y_sinh,y_cos,y_sin,y_dg,name,save]),
                 VBox([b_f,bbar_f,bbar_fbar,b_fbar,u_f,u_fbar,sum_f,sum_fb]),
                ])
amp_pars = interactive_output(plot_amplitudes,{
                'dm':dm,'dg':dg,'gs':gs,'r':r,'delta':delta,'gamma':gamma,'beta':beta,'afs':afs,'rez':rez,'imz':imz,
                'omega':omega,'d_omega':d_omega,'eff':eff,'d_eff':d_eff,'a_prod':a_prod,'a_det':a_det,'sigma_t':sigma_t,
                'xmin':xmin,'xmax':xmax,'y_cosh':y_cosh,'y_sinh':y_sinh,'y_cos':y_cos,'y_sin':y_sin,'y_dg':y_dg,
                'a_acc':a_acc,'n_acc':n_acc,'b_acc':b_acc,'beta_acc':beta_acc,'cutoff_acc':cutoff_acc,'k_acc':k_acc,
                'name':name,'save':save,
                'b_f':b_f,
                'bbar_f':bbar_f,
                'bbar_fbar':bbar_fbar,
                'b_fbar':b_fbar,
                'u_f':u_f,
                'u_fbar':u_fbar,
                'tagging':tagging,
                'k_tag':k_tag,
                'k_sum_f':sum_f,
                'k_sum_fb':sum_fb})
