#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file               :docker
@Description:       : python调用docker相关类
@Date               :2023/1/5 19:19:02
@Author             :lyzeng
@mail               :pylyzeng@gmail.com       
@version            :1.0
'''
import docker
from docker.tls import TLSConfig
from pathlib import Path
from dataclasses import dataclass, field

'''
docker -H=tcp://remote-server-ip:2376 --tlsverify --tlscacert=ca.pem --tlscert=cert.pem --tlskey=key.pem version
docker daemon --tls --tlscert=/path/to/server-cert.pem --t
'''

'''
# 语法参考使用：
mypassphrase = 'TLSdocker'
client = docker.DockerClient(base_url='tcp://remote-server-ip:2376', tls=True,
                            tls_client_key='ca-key.pem', tls_client_key_passphrase=mypassphrase,
                            tls_client_cert='ca.pem')
# 配置TLS连接
tls_config = TLSConfig(
    client_cert=('/path/to/client-cert.pem', '/path/to/client-key.pem'),
    ca_cert='/path/to/ca.pem',
    verify=True
)
# 连接到Docker守护进程
client = docker.DockerClient(base_url='tcp://localhost:2376', tls=tls_config, tls_client_key_passphrase=mypassphrase)
# 执行docker操作
project_name = "project_name"
number_of_solutions = "number_of_solutions"
best_PDB = "best_PDB"
restraint = "restraint"

output = client.containers.run(
    "pegi3s/pydock3",
    "bash -c './run_ftdock {} {} {} {}'".format(project_name, number_of_solutions, best_PDB, restraint),
    volumes={'/your/data/dir': {'bind': '/data', 'mode': 'rw'}},
    remove=True,
)
print(output)
'''
@dataclass()
class DockerClient():
    '''
    docker client
    # 配置TLS连接docker client
    tls_config = TLSConfig(
        client_cert=('/path/to/client-cert.pem', '/path/to/client-key.pem'),
        ca_cert='/path/to/ca.pem',
        verify=True
    )
    '''
    base_url: str # for mmme server is 'tcp://192.168.100.110:12376'
    # tls: TLSConfig = field(default=None)
    # passphrase: str = field(default=None)

    def __post_init__(self):
        # self.client = docker.DockerClient(base_url=self.base_url, tls=self.tls, tls_client_key_passphrase=self.passphrase)
        self.client = docker.DockerClient(base_url=self.base_url)

    def run(self, image="pegi3s/pydock3", command, volumes, remove=True):
        '''
        run docker container and return output(rm container after run)
output = client.containers.run(
    "pegi3s/pydock3",
    "bash -c './run_ftdock {} {} {} {}'".format(project_name, number_of_solutions, best_PDB, restraint),
    volumes={'/your/data/dir': {'bind': '/data', 'mode': 'rw'}},
    remove=True,
)
print(output)
        '''
        output = self.client.containers.run(
            image,
            command,
            volumes=volumes,
            remove=remove
        )
        return output
