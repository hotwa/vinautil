#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :comm
@Description:       :
@Date               :2023/1/2 21:26:21
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''
from json import dumps
import pickle
from pathlib import Path
import os

def piksave(obj:object,filename:str):
    with open(filename, 'wb') as f:
        pickle.dump(obj,f)

def pikload(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)
def bejson(d:dict(help='need to beauty json')) -> dict:
    return dumps(d,indent=4,ensure_ascii=False)



def find_file(env_var, keyword):
    # 获取环境变量中的所有路径
    paths = os.environ.get(env_var).split(os.pathsep)
    # 存储匹配的文件
    matches = []
    # 遍历每个路径
    for path in paths:
    # 判断路径是否存在
        if Path(path).exists():
            # 使用glob或rglob方法匹配包含关键字的文件，并将结果添加到匹配列表中（注意转义特殊字符）
            matches.extend(Path(path).glob(f"*{keyword}*"))
    # 判断匹配列表的长度
    if len(matches) == 0: # 没有找到任何文件，返回None
        return None
    elif len(matches) == 1: # 只找到一个文件，返回Path对象（取第一个元素）
        return matches[0]
    else: # 找到多个文件，抛出异常，并显示匹配列表中的所有文件名和路径信息
        raise Exception(f"Multiple files found: {matches}")