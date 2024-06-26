<script setup>
import axios from 'axios';
import { useNode, Handle, Position, useVueFlow, useHandleConnections } from '@vue-flow/core';
import { computed, watch } from 'vue';
const { node, connectedEdges } = useNode();
const { findNode, addEdges } = useVueFlow();
node.active = computed(() => {
    for (let edge of connectedEdges.value) {
        if (edge.target === node.id && !edge.animated) {
            return false;
        }
    }
    return true;
});
// useHandleConnections({
//   type: "source",
//   onConnect: (params) => {
//     params['id'] = `${params.source}-${params.target}`;
//     addEdges(params);
//     const source = findNode(params.source)
//     const edge = findEdge(params.id)
//     edge.animated = source.active;
//   },
// });

watch(() => node.active,
    (active) => {
        for (let edge of connectedEdges.value) {
            if (edge.source === node.id) {
                edge.animated = active;
            }
        }
    },
    {
        deep: true,
        immediate: true,
    }
)

async function runJob() {
    try{
        for (let key in node.data.taskIds) {
            const taskId = node.data.taskIds[key];
            if (taskId == null) {
                continue;
            }
            const response = await axios.get(base_url + "/inspect/" + taskId);
            if (response.data != "SUCCESS") {
                alert("job: " + taskId + " is " + response.data);
                return;
            }
        }
        const response = await axios.putForm(base_url + "/runJob/" + node.id, node.data.values);
        for (let rdt of response.data) {
            let found = false;
            for (let edge of connectedEdges.value) {
                if (edge.target === node.id) {
                    continue;
                }
                const targetNode = findNode(edge.target);
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

const base_url = import.meta.env.BASE_URL;
</script>

<template>
    <div>
        <span>{{ node.id }} </span>
        <button @click="runJob" :disabled="!node.active">run</button>
        <Handle type="source" :position="Position.Right" :is-valid-connection="(connection) => node.to.includes(connection.target)" />
        <Handle type="target" :position="Position.Left" :is-valid-connection="(connection) => node.from.includes(connection.source)" />
    </div>
</template>