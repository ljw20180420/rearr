<script setup>
import axios from 'axios';
import { useNode, Handle, Position, useVueFlow, useHandleConnections } from '@vue-flow/core';
import { ref, watch } from 'vue';
const { node } = useNode();
const { findNode, findEdge } = useVueFlow();
const connections = useHandleConnections({type: 'source'});
let targetEdges = [];
let targetNodes = [];
for (let i = 0; i < connections.value.length; ++i) {
    targetEdges.push(findEdge(connections.value[i].edgeId));
    targetNodes.push(findNode(targetEdges[i].target));
}
const sourceConns = useHandleConnections({type: 'target'})
let sourceEdges = [];
for (let sourceConn of sourceConns.value) {
    sourceEdges.push(findEdge(sourceConn.edgeId));
}

const disabled = ref(null);

watch(
    () => {
        for (let sourceEdge of sourceEdges) {
            if (!sourceEdge.animated) {
                return false;
            }
        }
        return true;
    },
    (activate) => {
        for (let targetEdge of targetEdges) {
            targetEdge.animated = activate;
        }
        disabled.value = !activate;
    },
    {
        deep: true,
        immediate: true
    }
)

async function runJob() {
    try{
        for (let key in node.data.taskIds) {
            const taskId = node.data.taskIds[key];
            if (taskId == null) {
                continue;
            }
            const response = await axios.get("inspect/" + taskId);
            if (response.data != "SUCCESS") {
                alert("job: " + taskId + " is " + response.data);
                return;
            }
        }
        const response = await axios.putForm("/runJob/" + node.id, node.data.values);
        for (let rdt of response.data) {
            let found = false;
            for (let targetNode of targetNodes) {
                for (let obj of targetNode.data) {
                    if (rdt.name == obj.name) {
                        obj.value = rdt.value;
                        obj.taskId = rdt.taskId;
                        found = true;
                    }
                    if (found) {
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
        }
    } catch(error) {
        alert(error);
    }
}
</script>

<template>
    <div>
        <span>{{ node.id }} </span>
        <button @click="runJob" :disabled="disabled">run</button>
        <Handle type="source" :position="Position.Right" />
        <Handle type="target" :position="Position.Left" />
    </div>
</template>