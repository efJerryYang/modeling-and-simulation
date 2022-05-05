1. `initialize!()`内部

   `gR[master_node]=rand(1:NUMNODE)`是可选的随机事件，代码中写为`gR[master_node]=1`，因为初始时是否随机选择主节点理论上不影响结果

2. `estimate_switch_state!()`内部

   - `rand!(distrX, lifeX)`，其中`X=A,B`，这里是从指数分布中取样，获得各个开关的寿命，结果写入`lifeX`向量中

   - `tolX=rand() * (1 - PX0)`, 其中`X=A,B`，这里是在开关已经损坏的条件下，随机选择一种损坏方式

3. `estimate_node_role_state!()`内部

   `alert_count = rand(NUM_NODE)`，具体为获取各个节点在失去主节点信号的时候警戒时间长度，等待时间最短的节点将被选为主节点

目前看来，好像就这三个，但在选节点的角色状态的时候可能还要有，但我这里本来就没有分析清楚。