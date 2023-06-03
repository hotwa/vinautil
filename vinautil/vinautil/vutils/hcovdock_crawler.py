#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :hcovdock_crawler
@Description:       : 爬虫分子对接HcovDock [文献](https://pubmed.ncbi.nlm.nih.gov/36573474/) [10.1093/bib/bbac559](https://doi.org/10.1093/bib/bbac559)
@Date               :2023/2/23 10:21:53
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''

import asyncio
from playwright.async_api import async_playwright
from playwright.sync_api import Playwright, sync_playwright, expect
from pathlib import Path
from time import sleep, strftime
import pickle



def submit_one_task(ligand_file: Path, receptor_file: Path, residue: str,
                    receptor_mail:str='pylyzeng@gmail.com',url="http://huanglab.phys.hust.edu.cn/hcovdock", header_string = "Task"):
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=False)
        context = browser.new_context()
        page = context.new_page()
        page.goto(url, timeout=100000)
        page.wait_for_load_state('networkidle')

        # 上传receptor
        input_locator = page.locator('#receptor')
        input_locator.set_input_files(receptor_file.as_posix())

        # 等待上传完成
        page.wait_for_selector('#residue_option')
        sleep(3)

        # 定位下拉框元素
        select_locator = page.select_option('#residue_option', residue)
        print(select_locator)

        # 上传ligand
        input_locator = page.locator('#ligand')
        input_locator.set_input_files(ligand_file.as_posix())

        # 等待上传完成
        page.wait_for_selector('#residue_option')
        sleep(3)

        page.fill('#emailaddress', receptor_mail)

        date_str = header_string + strftime("%Y%m%d")
        taskname = date_str + receptor_file.stem.split('.')[0]
        page.locator("input[name=\"jobname\"]").click()
        page.locator("input[name=\"jobname\"]").fill(taskname)
        page.get_by_role("checkbox").check()
        page.get_by_role("button", name="Submit").click()
        print(
            f''' Task submit:
taskname: {taskname}
receptor: {receptor_file.name}
ligand: {ligand_file.name}
mail: {receptor_mail} '''
        )




if __name__ == '__main__':
    # asyncio.run(main())
    root_path = Path('C:\\Users\\Admin\\OneDrive\\2020级研曾令宇\\实验数据\\hcovdock_benchmark')
    l_path = root_path.joinpath('molecular_polar_pdb')
    r_path = root_path.joinpath('receptor_pdb_scar')
    lfiles = list(l_path.glob('*.pdb'))
    rfiles = list(r_path.glob('*.scar.pdb'))
    # 读取之前提交的任务信息
    record_file = Path('task_record.pkl')
    if record_file.exists():
        with open(record_file, 'rb') as f:
            submitted_tasks = pickle.load(f)
    else:
        submitted_tasks = []
    for r,l in zip(rfiles, lfiles):
        t = r.stem.split('.')[0].split('-')[1].split('_')[1:3]
        target_residue = f"{t[0]}:GLY{t[1]}"
        task_info = (r.name, l.name, target_residue)
        if task_info in submitted_tasks:
            print(f"Task {task_info} has already been submitted, skipped.")
            continue
        submit_one_task(ligand_file=l, receptor_file=r, residue=target_residue, receptor_mail='pylyzeng@gmail.com', header_string='Newtask_gmail-')


        # 将新提交的任务信息记录下来
        submitted_tasks.append(task_info)
        with open(record_file, 'wb') as f:
            pickle.dump(submitted_tasks, f)

        sleep(8)
