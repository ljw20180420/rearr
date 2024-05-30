<script setup>
import { ref, watch } from 'vue';
import axios from 'axios';
import { useNode, Handle, Position, useVueFlow, useHandleConnections } from '@vue-flow/core';
const files = ref(null);
const percentage = ref(0);
const { node, connectedEdges } = useNode();
const { findNode, addEdges } = useVueFlow();
useHandleConnections({
    type: "source",
    // onConnect: (params) => {
    //     for (let param of params) {
    //         param['id'] = `${param.source}-${param.target}`;
    //         addEdges(param);
    //         const source = findNode(param.source)
    //         const edge = findEdge(param.id)
    //         const target = findNode(param.target)
    //         let activate = true;
    //         for (let obj of source.data) {
    //             target.data.values[obj.name] = obj.value
    //             if (Object.keys(obj).includes("taskId")) {
    //                 target.data.taskIds[obj.name] = obj.taskId
    //             }
    //             if (obj.value == null || obj.value === "" ) {
    //                 activate = false;
    //             } else if (Array.isArray(obj.value)) {
    //                 for (let val of obj.value) {
    //                     if (val == null) {
    //                         activate = false;
    //                     }
    //                 }
    //             }
    //             edge.animated = activate;
    //         }
    //     }
    // },
    onDisconnect: (params) => {
        console.log("onDisConnect", params);
        for (let param of params) {
            const source = findNode(param.source)
            const target = findNode(param.target)
            for (let obj of source.data) {
                if (Array.isArray(obj.value))
                {
                    target.data.values[obj.name] = Array(obj.value.length).fill(null);
                } else {
                    target.data.values[obj.name] = null;
                }
                if (Object.keys(obj).includes("taskId")) {
                    target.data.taskIds[obj.name] = null
                }
            }
        }
    },
});

async function upload(event, idx) {
    files.value = event.target.files;
    try{
        const response = await axios.putForm("/upload", {
            files: files.value,
        },
        {
            onUploadProgress : progressEvent => {
                percentage.value = Math.round((progressEvent.loaded * 100) / progressEvent.total);
            }
        });
        node.data[idx].taskId = response.data.taskId;
        for (let i in node.data[idx].value) {
            node.data[idx].value[i] = response.data.value[i];
        }
    } catch(error) {
        alert(error);
    }
}

watch(() => node.data,
    (data) => {
        for (let edge of connectedEdges.value) {
            if (edge.target === node.id) {
                continue;
            }
            const targetNode = findNode(edge.target);
            let activate = true;
            for (let obj of data) {
                targetNode.data.values[obj.name] = obj.value
                if (Object.keys(obj).includes("taskId")) {
                    targetNode.data.taskIds[obj.name] = obj.taskId
                }
                if (obj.value == null || obj.value === "") {
                    activate = false;
                } else if (Array.isArray(obj.value)) {
                    for (let val of obj.value) {
                        if (val == null) {
                            activate = false;
                        }
                    }
                }
            }
            edge.animated = activate;
        }
    },
    {
        deep: true,
        immediate: true,
    }
)
</script>

<template>
    <div>
        <div v-for="obj, i in node.data">
            <span>{{ obj.name }}: </span>
            <a v-if="Object.keys(obj).includes('taskId')" :href="'/download/' + obj.taskId">{{ obj.taskId }}</a>
            <br>
            <input v-if="obj.type == 'value'" v-model="node.data[i].value">
            <div v-else-if="obj.type == 'file' || obj.type == 'files'">
                <input type="file" @change="upload($event, i)" :multiple="obj.type == 'files'">
                <span>{{ percentage }}</span>
            </div>
            <select v-else-if="obj.type == 'select'" v-model="node.data[i].value">
                <option v-for="option of node.data[i].options">
                    {{ option }}
                </option>
            </select>
            <hr>
        </div>
        <Handle type="source" :position="Position.Right" :is-valid-connection="(connection) => node.to.includes(connection.target)" />
        <Handle type="target" :position="Position.Left" :is-valid-connection="(connection) => node.from.includes(connection.source)" />
    </div>
</template>