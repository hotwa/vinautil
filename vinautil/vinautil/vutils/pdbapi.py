#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :pdbapi
@Description:       :
@Date               :2023/6/26 18:25:07
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''
from dataclasses import dataclass, field


@dataclass
class PDBAPI:
    """https://data.rcsb.org/redoc/
    """
    # entry_id: str = field(default='1a4w', metadata={'help': 'PDB entry id'})
    # asym_id: str = field(default='A', metadata={'help': 'asym id'})
    comp_id: str = field(default='ZYA', metadata={'help': 'comp id'})

    def querymole():
        # 定义API的URL
        url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}"
        # 发送GET请求
        response = requests.get(url)
        # 将返回的JSON数据转换为Python字典
        data = response.json()
        return data