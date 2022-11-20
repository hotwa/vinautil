#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file        :pymol_color_plugin.py
@Description:       :pymol配色方案
@Date     :2021/11/10 10:52:25
@Author      :hotwa
@version      :1.0
'''


from pymol import cmd
from pathlib import Path

class color_plan():
    color1 = (('color1', '[186,182,217]'), 'purple')
    color2 = (('color2', '[233,195,153]'), 'yellow')
    color3 = (('color3', '[43,113,216]'), 'blue-N')
    color4 = (('color4', '[206,155,198]'), 'purple')
    color5 = (('color5', '[251,187,62]'), 'orange')
    color6 = (('color6', '[245,157,158]'), 'red')
    color7 = (('color7', '[133,188,135]'), 'green')
    color8 = (('color8', '[30,230,30]'),'green-CL') # Cl卤素配色
    color9 = (('color9', '[141,215,247]'),'blue-C') # C配色
    color10 = (('color10', '[0,132,55]'),'green-F') # F卤素配色
    grey1 = ('grey1','[224,224,224]')
    colors = (color1, color2, color3, color4, color5, color6, color7, color8, color9, color10)

    def __init__(self, path, id):
        self.id = id.lower()
        cmd.reinitialize()
        p = Path(path)
        if p.is_dir():
            file = p.joinpath(f"{id}.pdb")
        else:
            raise ValueError('path params error')
        cmd.load(file,id)

    def pretty(self):
        cmd.remove('solvent metals') # 移除金属离子和溶剂
        cmd.pretty(selection=self.id)

    @staticmethod
    def defcolor():
        # 定义常用配色
        color_gennerate = map(lambda x:x[0],color_plan.colors)
        list(map(lambda x:cmd.set_color(*x),color_gennerate))
        # 定义灰度配色
        cmd.set_color(*color_plan.grey1)

    @staticmethod
    def grey(obj='sele',opt=True):
        color_plan.defcolor()
        if opt: method.optimisation() # 优化
        # 对某个对象进行灰化
        cmd.color('grey1',f'{obj}')

    @staticmethod
    def color_mole(obj='hetatm'):
        # color blue,sele in name c*
        color_plan.defcolor()
        cmd.color('color9',f'{obj} in name c*')
        cmd.color('color6',f'{obj} in name o*')
        cmd.color('color5',f'{obj} in name s*')
        cmd.color('color3',f'{obj} in name n*')
        cmd.color('color8',f'{obj} in name cl*')
        cmd.color('color10',f'{obj} in name f*')
    
    @staticmethod
    def pretty(obj = None):
        if not obj: obj = cmd.get_names()[0] 
        cmd.remove('solvent metals') # 移除金属离子和溶剂
        cmd.pretty(selection=obj)
        

class font():
    font_size = 28 # 单位就是正常的px。你也可以用负值，则单位是Å

class mole():
    stick_radius=0.10

class atom():
    size = 0.28

class method():

    @staticmethod
    def optimisation(grey=True):
        color_plan.defcolor()
        print('常规优化')# 常规优化
        if grey: cmd.set('cartoon_color','grey1')
        cmd.set('stick_transparency','0.1')
        cmd.set("ray_shadows","off") # 消除投射阴影
        cmd.set('cartoon_highlight_color') # 增加不可见区域
        cmd.set('cartoon_fancy_helices') # β折叠优化
        cmd.set('cartoon_transparency','0.5') # 透明度设置0.5

    @staticmethod
    def remove():
        cmd.remove('solvent metals') # 移除金属离子和溶剂
        

cmd.extend('color_plan',color_plan)
cmd.extend('method',method)