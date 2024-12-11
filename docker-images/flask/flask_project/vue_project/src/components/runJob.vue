<script setup>
import axios from 'axios';
import { useNode, Handle, Position, useVueFlow } from '@vue-flow/core';
import { computed } from 'vue';
import { isValidConnection } from './utils.js'
const { node, connectedEdges } = useNode();
const { findNode } = useVueFlow();

node.active = computed(() => {
    for (let edge of connectedEdges.value) {
        if (edge.target === node.id && !edge.animated) {
            return false;
        }
    }
    return true;
});

node.data = computed(() => {
    let data = { values: {}, taskIds: {} };
    for (let edge of connectedEdges.value) {
        if (edge.source === node.id) {
            continue;
        }
        for (let obj of edge.sourceNode.data) {
            data.values[obj.name] = obj.value
            if (Object.keys(obj).includes("taskId")) {
                data.taskIds[obj.name] = obj.taskId;
            }
        }
    }
    return data;
})

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
    <div :title="node.title">
        <span>{{ node.id }} </span>
        <button @click="runJob" :disabled="!node.active">run</button>
        <Handle id="s" type="source" :position="Position.Right" :is-valid-connection="isValidConnection" />
        <Handle id="t" type="target" :position="Position.Left" :is-valid-connection="isValidConnection" />
    </div>
</template>